/**
 * Project iDynoMicS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * This solver uses a fast fourier transform (fft) algorithm to solve the
 * diffusion equation explicitly. It is unconditionally stable, but has some
 * restriction that should be remembered: - grid size has to be an integer power
 * of 2 (2^n, e.g. 32,64,128). other algorithms can be implemented - fft is
 * intrinsically periodic, hence non-periodic boundaries are not (yet, TODO)
 * supported. - reaction rates can be included, if they are assumed to be
 * (quasi-)constant at the used timescale
 * @since June 2006
 * @version 0.1
 * @author Andreas D�tsch, 
 */

package simulator.diffusionSolver;

import idyno.SimTimer;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.SpatialGrid;
import simulator.geometry.DiscreteVector;
import simulator.reaction.Reaction;
import utils.FftOperations;
import utils.LogFile;
import utils.XMLParser;

/**
 * TODO: implement diffusion field (inhomogeneous diffusion)
 */
public class Solver_FFT extends DiffusionSolver {

	/* +++++++++++ NEEDS TESTING +++++++++++ */
	private static final long     serialVersionUID = 1L;

	// size of the original grid
	protected int                 _nI, _nJ, _nK;
	protected int                 _ntot;
	protected int[]               _gridSize        = new int[3];
	protected double              _lX, _lY, _lZ;
	public SoluteGrid			  bioMass;
	// grid resolution
	protected double              _dx;
	// time step
	protected double              _dT;

	// Grids storing the current reaction rates of the solutes
	// protected ArrayList<SoluteGrid> _rates = new ArrayList<SoluteGrid>();
	protected SoluteGrid[]        _rateGrids,_diffRateGrids;

	// Padding is not used for FFT, but could be to implement non-periodic
	// geometries (TODO)
	protected boolean             _usesPadding     = false;

	// --- auxiliary variables ---

	// wave number sum of squares
	protected double[][][]        _kSquareSum;

	// temporary variables for storage of grids as 1D array
	protected double[]            _vectorizedSolute;
	protected double[]            _vectorizedRate;

	// TODO remove the DiscreteCoordinate once the getValueAt methods are
	// corrected for padding-handling
	protected DiscreteVector      _currentLocation = new DiscreteVector();

	protected boolean             _isInitialized   = false;

	// epsilon, the smallest number that is non-zero
	protected static final double _E               = 2.220446049250313e-016;
	private static int            nSolute, currentSolute;

	/**
     * 
     */
	public void init(Simulator aSimulator, XMLParser aSolverRoot, Double agentTimeStep) {
		super.init(aSimulator, aSolverRoot, agentTimeStep);
		initializeConcentrationFields();
		nSolute = _soluteList.length;
	}

	/**
     * 
     */
	public void initializeConcentrationFields() {
		// Re-initialization is not necessary since the grid doesn't change
		if (_isInitialized) { return; }

		// TODO make sure that all solute grids are of the same size

		// initialize grid size
		_nI = _soluteList[0].getGridSizeI();
		_nJ = _soluteList[0].getGridSizeJ();
		_nK = _soluteList[0].getGridSizeK();
		_gridSize[0] = _nI;
		_gridSize[1] = _nJ;
		_gridSize[2] = _nK;
		_ntot = _nI*_nJ*_nK;
		_dx = _soluteList[0].getResolution();
		_lX = _nI*_dx;
		_lY = _nJ*_dx;
		_lZ = _nK*_dx;
		bioMass = new SoluteGrid(_nI, _nJ, _nK, _dx);
		
		
		// initialize temporary storage grid for reaction rates
		_rateGrids = new SoluteGrid[_soluteList.length];
		_diffRateGrids = new SoluteGrid[_soluteList.length];
		for (int i = 0; i<_rateGrids.length; i++)
		{
			_rateGrids[i] = new SoluteGrid(_nI, _nJ, _nK, _dx);
			_diffRateGrids[i] = new SoluteGrid(_nI, _nJ, _nK, _dx);
		}
		

		// initialize temporary vectorized grid. FFT works with 1D arrays of
		// complex numbers, so n-D arrays of real numbers are stored as vectors
		// with alternating real and imaginary parts.
		_vectorizedSolute = new double[2*_nI*_nJ*_nK];
		_vectorizedRate = new double[2*_nI*_nJ*_nK];

		// initialize wave number representation for fft
		double[] kX = new double[_nI];
		double[] kY = new double[_nJ];
		double[] kZ = new double[_nK];
		double qX = 2*Math.PI/_lX;
		double qY = 2*Math.PI/_lY;
		double qZ = 2*Math.PI/_lZ;
		for (int i = 0; i<_nI; i++) {
			kX[i] = (i<=_nI/2 ? qX*i : qX*(i-_nI));
		}
		for (int j = 0; j<_nJ; j++) {
			kY[j] = (j<=_nJ/2 ? qY*j : qY*(j-_nJ));
		}
		for (int k = 0; k<_nK; k++) {
			kZ[k] = (k<=_nK/2 ? qZ*k : qZ*(k-_nK));
		}

		_kSquareSum = new double[_nI][_nJ][_nK];
		for (int i = 0; i<_nI; i++) {
			for (int j = 0; j<_nJ; j++) {
				for (int k = 0; k<_nK; k++) {
					_kSquareSum[i][j][k] = kX[i]*kX[i]+kY[j]*kY[j]+kZ[k]*kZ[k];
				}
			}
		}
		// set (0,0,0) to epsilon, to prevent division by zero error
		_kSquareSum[0][0][0] = _E;

		// check system validity. The used FFT only works with integer power of
		// 2.
		try {
			validateSystem();
		} catch (Exception e) {
			System.out.println("Error: All dimensions have to be an integer power of 2!");
			System.exit(1);
		}
		_isInitialized = true;
	}

	/**
     * 
     */
	public void solveDiffusionReaction() {
		
		_dT = SimTimer.getCurrentTimeStep();

		// update reaction rates
		computeReactionRates();

		// solve for all solutes

		for (int iChem = 0; iChem<nSolute; iChem++) {
			
			LogFile.writeLogAlways("Start Solute "+iChem);
			
			currentSolute = /*_soluteIndex.get(*/ iChem /*)*/ ;
			LogFile.writeLogAlways("currentSolute "+currentSolute);
			// make a 1-D alternating grid
			vectorizeGrid(currentSolute);			
			// transform grids, store in temporary vectorized arrays
			_vectorizedSolute = FftOperations.fftn(_vectorizedSolute, 3, _gridSize, 1);			
			_vectorizedRate = FftOperations.fftn(_vectorizedRate, 3, _gridSize, 1);			
			// Calculate solution
			solveInFrequencySpace(currentSolute);			
			// inverse transformation (normalization is postponed)
			_vectorizedSolute = FftOperations.fftn(_vectorizedSolute, 3, _gridSize, -1);
			// rewrite 3D grid
			updateSoluteGrid(currentSolute);
			
			LogFile.writeLogAlways("Solute "+iChem+": done");
			
		}
		LogFile.writeLogAlways("Solutes finished");
	}

	/**
     * Transforms a 3D spatialGrid to a 1D "vectorizedGrid", suitable for
     * fft-algorithm. Values of the original Grid are stored by rows, i.e. the
     * first subscript (i) changes least rapidly. Each value occupies 2 values
     * in the vectorized grid, representing the real and imaginary part of a
     * complex number. Since the input is real by default, the imaginary parts
     * are initialized with 0.
     * 
     * @param iSolute index of the soluteGrid to transform
     */
	protected void vectorizeGrid(int iSolute) {
		int ii = 0;
		for (int i = 0; i<_nI; i++) {
			for (int j = 0; j<_nJ; j++) {
				for (int k = 0; k<_nK; k++) {
					// get solute grid WITHOUT PADDING
					_currentLocation.set(i, j, k);
					_vectorizedSolute[ii] = _soluteList[iSolute].getValueAt(_currentLocation);
					_vectorizedSolute[ii+1] = 0.0D;
					_vectorizedRate[ii] = _rateGrids[iSolute].getValueAt(_currentLocation);
					_vectorizedRate[ii+1] = 0.0D;
					ii += 2;
				}
			}
		}
	}

	/**
     * Inversion of vectorizeGrid. Rewrites the stored-in-rows values of a
     * vectorized grid to the original 3D array.
     * 
     * @param iSolute
     */
	protected void updateSoluteGrid(int iSolute) {
		int ii = 0;
		for (int i = 0; i<_nI; i++) {
			for (int j = 0; j<_nJ; j++) {
				for (int k = 0; k<_nK; k++) {
					_currentLocation.set(i, j, k);
					_soluteList[iSolute].setValueAt(_vectorizedSolute[ii]/_ntot, _currentLocation);
					ii += 2;
				}
			}
		}
	}

	protected void solveInFrequencySpace(int iSolute) {

		// auxiliary temp variables
		int ii = 0;
		double D = _soluteList[iSolute].getDiffusivity();
		double[] zeroK = new double[2];
		double tempDK = 0;
		double tempExp = 0;

		// solve for k=(0,0,0)
		zeroK[0] = _vectorizedRate[0]*_dT+_vectorizedSolute[0];
		zeroK[1] = _vectorizedRate[1]*_dT+_vectorizedSolute[1];

		// solve for any point in frequency space
		for (int i = 0; i<_nI; i++) {
			for (int j = 0; j<_nJ; j++) {
				for (int k = 0; k<_nK; k++) {
					tempDK = D*_kSquareSum[i][j][k];
					tempExp = Math.exp(-_dT*tempDK);

					_vectorizedSolute[ii] = _vectorizedSolute[ii]*tempExp+_vectorizedRate[ii]
					        *(1-tempExp)/tempDK;
					_vectorizedSolute[ii+1] = _vectorizedSolute[ii+1]*tempExp+_vectorizedRate[ii+1]
					        *(1-tempExp)/tempDK;

					ii += 2;
				}
			}
		}

		// remember pre-computed value for k=(0,0,0)
		_vectorizedSolute[0] = zeroK[0];
		_vectorizedSolute[1] = zeroK[1];
	}

	/**
     * update reaction rate fields for each solute according to the current
     * reaction rates. reaction rates for solutes are called from the reactions
     * and stored as a reaction grid here
     */
	protected void computeReactionRates() {

		// first reset the grid
		for (int i = 0; i<_rateGrids.length; i++) {
			_rateGrids[i].setAllValueAt(0);
		}

		// Now add on each grid the contribution of each reaction
		
		LogFile.writeLogAlways("initializing counters");
		
		for (int i=0; i<_nI; i++)
		{	for (int j=0; j<_nJ; j++)
			{	for (int k=0; k<_nK; k++)
				{
					LogFile.writeLogAlways("counters set:"+i+" "+j+" "+k);
					// FFT has no biomass-Grid -> created one for a walkaround
					bioMass.grid[i][j][k]=0.0;
				}
			}
		}
					LogFile.writeLogAlways("bioMass resetted");
					for (Reaction aReaction : _reactions) 
					{
						for (int p=0; p<_reactions.size(); p++)
							{ _reactions.get(p).fitAgentMassOnGrid(bioMass); }
						aReaction.applySignalReaction(_soluteList, _rateGrids,_diffRateGrids, bioMass);
					}
	}

	/**
     * Check the system validity. The fft algorithm can only handle arrays of a
     * length that is an integer power of 2. If one of the dimensions does not
     * fit, an exception is thrown.
     * 
     * @throws NoPowerOfTwoException
     */
	public void validateSystem() throws Exception {
		boolean allPowerOfTwo = true;
		allPowerOfTwo &= utils.FftOperations.checkPowerOfTwo(_nI);
		allPowerOfTwo &= utils.FftOperations.checkPowerOfTwo(_nJ);
		allPowerOfTwo &= utils.FftOperations.checkPowerOfTwo(_nK);
		if (!allPowerOfTwo) throw new Exception();
	}
}
