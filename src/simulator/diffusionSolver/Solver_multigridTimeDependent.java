
/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 *  
 *______________________________________________________
 * DiffusionSolver is an abstract class used as parent for all diffusion_solvers 
 * you could define
 * 
 * This class lets you solve a time dependent nonlinear reaction diffusion system.
 * Implemented is Crank-Nicolson in combination with Full Approximation Scheme (FAS, borrowed from class Solver_multigrid).
 * 
 */

/**
 * @since August 2011
 * @version 1.0
 * @author Katrin Bohl, FSU Jena
 * @author Max Adolph, TU Dresden
 */
package simulator.diffusionSolver;

import java.util.ArrayList;
import idyno.SimTimer;
import simulator.diffusionSolver.multigrid.MultigridSoluteTimeDependent;
import simulator.geometry.Domain;
import simulator.geometry.boundaryConditions.AllBC;
import simulator.Simulator;
import simulator.SoluteGrid;
import utils.LogFile;
import utils.XMLParser;

public class Solver_multigridTimeDependent extends DiffusionSolver
{
	
	protected MultigridSoluteTimeDependent 		_bLayer, _diffusivity;
	
	protected MultigridSoluteTimeDependent[] 	_solute, _biomass;
	
	protected SoluteGrid[]      				allSolute, allReac, allDiffReac; 
	
	/**
	 * Stores the reaction rates for one single reaction (without the
	 * catalyst/biomass).
	 */
	public SoluteGrid[] 						singleReac;
	
	/**
	 * Stores the reactions of each reaction during the time-step (without the
	 * catalyst/biomass).
	 */
	public SoluteGrid[]							reacSum;
	
	/**
	 * The default solver is used by this solver to out-source fast solutes.
	 */
	public Solver_multigrid						steadyStateSolver;
	
	/**
	 * The agent-timestep.
	 */
	public Double 								timeStep;
	
	/**
	 * the smaller timestep, that determines the uptake-frequency TODO make it easily adjustable
	 */
	public Double 								uptakeStep = 0.01;

	/**
	 * Commonly used counter for nSolutes.
	 */
	protected static int        				iSolute;
	
	/**
	 * The number of solutes in the system.
	 */
	protected int								nSolute; 
	
	/**
	 * The number of reactions in our system.
	 */
	protected int								nReaction;
	
	protected Domain            				_domain; 
	
	/**
	 * Spatial resolution of our coarsest grid (might be the finest if
	 * Multigrid is not used).
	 */
	public Double h = 0.0;
	
	public void init(Simulator aSimulator, XMLParser xmlRoot, Double agentTimeStep) 
	{
		super.init(aSimulator, xmlRoot, agentTimeStep);
		
		timeStep=agentTimeStep;
		
		// Create the table of solute grids.
		nSolute = _soluteList.length;
		_solute = new MultigridSoluteTimeDependent[nSolute];
		singleReac = new SoluteGrid[_reactions.size()];	
		reacSum = new SoluteGrid[_reactions.size()];
		allSolute = new SoluteGrid[nSolute];
		allReac = new SoluteGrid[nSolute];
		allDiffReac = new SoluteGrid[nSolute];
		// Initial concentration for each solute solved by this solver.
		
		_domain = aSimulator.world.getDomain(xmlRoot.getAttribute("domain"));
		_bLayer = new MultigridSoluteTimeDependent(_soluteList[0], "boundary layer");
		_diffusivity = new MultigridSoluteTimeDependent(_soluteList[0], "relative diffusivity");
		LogFile.writeLogAlways("Made bLayer and diffusivity");
		for ( int i = 0; i < nSolute; i++ ) 
		{
			if ( _soluteIndex.contains(i) ) 
			{
				Double sBulk = mySim.world.getMaxBulkValue(_soluteList[i].soluteIndex);
				_solute[i] = new MultigridSoluteTimeDependent(_soluteList[i],
						_diffusivity, _bLayer, mySim, sBulk, _reactions.size());
			} 
			else 
			{
				_solute[i] = null;
			}
			LogFile.writeLogAlways("Sorted solute: "+i);
		}
		LogFile.writeLogAlways("AgentTimeStep:"+ aSimulator.agentTimeStep);
		/*
		 * From this moment, nSolute is the number of solutes SOLVED by THIS
		 * solver.
		 */
		nSolute = _soluteIndex.size();
		nReaction = _reactions.size();
		// Initialize array of reactive biomasses.
		_biomass = new MultigridSoluteTimeDependent[nReaction];
		for (int i = 0; i<nReaction; i++) 
		{
			_biomass[i] = new MultigridSoluteTimeDependent(_soluteList[0],
											_reactions.get(i).reactionName);
			_biomass[i].resetMultigridCopies(0.0);
			LogFile.writeLogAlways("Made biomass "+i);
		}
		LogFile.writeLogAlways("sorting solutes");
		/*
		 * TODO do a seperate Method
		 * computes the size of one voxel (real size / number of voxels)
		 */
		h = _solute[0]._referenceSystemSide;
		h /= _solute[0].referenceIndex(_solute[0]._conc.getGridSizeI(),
										_solute[0]._conc.getGridSizeJ(),
										_solute[0]._conc.getGridSizeK());
		/*
		 * Adjust _soluteIndex and create cutSolutes as the indices for
		 * steadyState.
		 */
		ArrayList<Integer> cutSolutes = sortSolutes();
		
		LogFile.writeLogAlways("creating second solver now");
		
		if ( !cutSolutes.isEmpty() ) //checks, if a solute had to be removed from this solver
		{
			LogFile.writeLogAlways("solutes cut");
			// all removed solutes are handled by a steady-state-solver
			try
			{
				steadyStateSolver = (Solver_multigrid) 
						Class.forName("simulator.diffusionSolver.Solver_multigrid").newInstance();
			}
			catch (Exception e)
			{
				LogFile.writeLogAlways("Failed to create Solver_multigrid: "+e);
			}
			// fills the empty solver
			steadyStateSolver.initCoupled(mySim, _domain, cutSolutes, _reactions);
			
			// shallow copy all solutes
			steadyStateSolver.nSolute=nSolute;
			LogFile.writeLogAlways("steadyStateSolver created");
		}
	}

	/**
	 * Removes all 'fast' solutes from the jurisdiction of this solver 
	 * and created a List of removed solutes.
	 *
	 * @return List of removed solutes (inverse of the resulting _soluteIndex).
	 */
	public ArrayList<Integer> sortSolutes()
	{
		
		ArrayList<Integer> cutSolutes = new ArrayList<Integer>();
		int indexShift = 0;
		LogFile.writeLogAlways("Solutes: "+nSolute);
		for ( int i = 0;i < nSolute; i++ )
		{
			Double D = _solute[i].realGrid.diffusivity;
			/*
			 * TODO sort out the padding-stuff in the Grids
			 * and then use getGridSizeI() etc.
			 */
			Double nX = (_solute[i]._conc.grid.length-2) * h;
			Double nY = (_solute[i]._conc.grid[1].length-2) * h;
			Double nZ = (_solute[i]._conc.grid[1][1].length-2) * h;
			
			LogFile.writeLogAlways("D: "+D+" / nX: "+nX+" / nY: "+nY+" / nZ: "+nZ);
			
			Double maxArea = Math.max(nX, nY);
			maxArea = Math.max(maxArea, nZ);
			maxArea = maxArea*maxArea;
			
			if ( (timeStep*D) >= maxArea ) //timeStep or uptakeStep?
			{
				// Terminate solute from soluteIndex.
				cutSolutes.add(_soluteIndex.remove(i - indexShift));
				indexShift++;
			}
			//LogFile.writeLogAlways("finished processing: "+i);
		}
		return cutSolutes;
	}
	
	public void initializeConcentrationFields() 
	{
		// double minimalTimeStep = SimTimer.getCurrentTimeStep()/10;

		// Refresh then insert here the boundary layer and the diffusivity grid
		_domain.refreshBioFilmGrids();

		_bLayer.setGrid(_domain.getBoundaryLayer());
		
		_diffusivity.setGrid(_domain.getDiffusivity());
		
		singleReac = new SoluteGrid[_reactions.size()];	
		reacSum = new SoluteGrid[_reactions.size()];
		
		for (int i=0;i<_reactions.size();i++)
		{
			singleReac[i] = new SoluteGrid(_bLayer._conc.grid.length, _bLayer._conc.grid[1].length, 
					_bLayer._conc.grid[1][1].length, _bLayer._conc._reso);
			reacSum[i] = new SoluteGrid(_bLayer._conc.grid.length, _bLayer._conc.grid[1].length, 
					_bLayer._conc.grid[1][1].length, _bLayer._conc._reso);
		}
		// Prepare a soluteGrid with catalyst CONCENTRATION
		for (int i = 0; i<_biomass.length; i++) 
			{ _biomass[i].resetFinest(0d);
			_reactions.get(i).fitAgentMassOnGrid(_biomass[i].getGrid()); }
		
		for (int iSolute : _soluteIndex)
			_solute[iSolute].readBulk();
	}
	

	/**
	 * Solve by iterative relaxation
	 */
	public void solveDiffusionReaction() 
	{
		double timeToSolve = SimTimer.getCurrentTimeStep();
		internalIteration = 0;
		internTimeStep = timeToSolve;

		// bvm note 13.7.09:
		// this iterative loop is only passed through once because of
		// the value of internTimeStep used above; we leave the loop
		// as-is though to allow future use of iterates if needed
		
		// Katrin 11.8.11: I used this for the time dependence loop
		while (timeToSolve>0) 
		{
			// solve the simulation of diffusion and reaction
			stepSolveDiffusionReactionNoMultiGrid();
			System.out.println("Diffusion calculation at iteration "+internalIteration);
			
			// update bulk concentration
			//updateBulk();

			// Manage iterations
			internalIteration += 1;
			timeToSolve = timeToSolve-internTimeStep;		
		}				
	}
	
	/**
	 * Solves our diffusion grid and applies the uptake in small timesteps to simulate our domain.
	 * Because we have solutes that are diffusing rapidly, we can assume them to reach
	 * steady state within one step of uptake, so we use the faster steady-state-solver for them.
	 * The other ones get solved timedependently. The uptake is done just for the time-dependent 
	 * soluted, because the steady-state-solver incorporates the reaction into its calculation
	 * by design. 
	 * 
	 * @author Max Adolph, TU Dresden
	 */
	public void coupledSolver() 
	{	
		LogFile.writeLogAlways("___start___");
		//the dimensions of the finest grid
		int xFine=_solute[0]._conc.grid.length;
		int yFine=_solute[0]._conc.grid[1].length;
		int zFine=_solute[0]._conc.grid[1][1].length;
		
		// set the uptakeStep back to default
		double CurrentUptakeStep = uptakeStep;
		
		double time=0.0; //counts the actual time that has been solved already
		boolean lastStep=false; //gets set to true when last step is reached and differs behaviour
		double soluteDiff; //the difference in solute made by uptake
		double fudgeFactor; //necessary alteration for small concentrations 
		
		// sets the dt for every solute independently
		for (int iSolute : _soluteIndex)
			_solute[iSolute].computeStableDt(CurrentUptakeStep);
		
		// execute the actual solver
		while (!lastStep) //handles the uptake-timescale
		{
			LogFile.writeLogAlways("time: "+time);
			if ((timeStep-time)<=(CurrentUptakeStep)) // in case it is the last step, the simulated time is scaled down
			{
				CurrentUptakeStep=timeStep-time;
				
				for (int iSolute : _soluteIndex) //adjust dt to the new uptakeStep (smaller than before)
					_solute[iSolute].computeStableDt(CurrentUptakeStep);
				
				lastStep=true;
			}
		
			// handles the diffusion-timescale
			for (int iSolute : _soluteIndex) // does one (uptake)-time-step for every solute in the domain
			{
				double solTime=0.0; //internal time for this solute
				while ((solTime+_solute[iSolute].simTimeStep)<=(CurrentUptakeStep))
				{
					solTime=solTime+_solute[iSolute].simTimeStep;
						
					LogFile.writeLogAlways("___tik-tok, what says the clock?: "+(solTime+time));
					//apply the implicit solver to the concentration-domain
					_solute[iSolute].implCrankNic(_solute[iSolute]._conc, _solute[iSolute].simTimeStep, 
								_bLayer, _diffusivity._conc.grid, true);
				}
				_solute[iSolute].applyComputation();
			}
			
			// solve fast solutes with the steady-state-solver
			if (steadyStateSolver != null) //check, if steady state solver actually was created
			{
				/**
				 * because we just quickly scripted the steady-state-Solver it did not initialise
				 * its own biomass, reldiff etc. However, he has to use our computed grids anyway,
				 * so we just cram our values into it.
				 */
				
				// set boundary-layer
				steadyStateSolver._bLayer._conc[steadyStateSolver.maxOrder-1].grid=_bLayer._conc.grid;
				steadyStateSolver._bLayer.restrictToCoarsest("bLayer");
				for (int i=0;i<nSolute;i++)
				{
					// set the conc-grids
					steadyStateSolver._solute[i]._conc[steadyStateSolver.maxOrder-1].grid=_solute[i]._conc.grid;
					steadyStateSolver._solute[i].realGrid.grid=_solute[i].realGrid.grid;
					// set the relative diffusivity
					steadyStateSolver._diffusivity._conc[steadyStateSolver.maxOrder-1].grid=_diffusivity._conc.grid;
					steadyStateSolver._diffusivity.restrictToCoarsest("relDiff");
					
					for (int p=0;p<_solute[i]._conc.grid.length;p++)
						LogFile.writeLogAlways("solute "+i+": "+_solute[i]._conc.grid[p][1][1]);
				}
				// set biomass-grids
				for (int i=0;i<nReaction;i++)
				{
					steadyStateSolver._biomass[i]._conc[steadyStateSolver.maxOrder-1].grid=_biomass[i]._conc.grid;
					steadyStateSolver._biomass[i].restrictToCoarsest("bLayer");
				}
				// run steady-state-solver
				steadyStateSolver.stepSolveDiffusionReactionPlugIn();
			}
			//our diffusion is now done
			
			// do uptake
			updateReacRateAndDiffRate(); //adjust parameters to new conditions
			for (int iSolute : _soluteIndex)
				for (int i=0;i<xFine;i++)
					for (int j=0;j<yFine;j++)
						for (int k=0;k<zFine;k++)
						{
							soluteDiff = _solute[iSolute]._reac.grid[i][j][k]*CurrentUptakeStep; //compute change in solute
							if ( -soluteDiff > _solute[iSolute]._conc.grid[i][j][k] ) //check for unreasonable uptake
							{
								LogFile.writeLogAlways("Solute "+iSolute+" is getting fudged");
								/*
								 * alterated computation of the uptake to avoid linearization-errors
								 * 
								 * the rate of consumption varies largely for small concentrations 
								 * and would rapidly deplete even during one uptake-step
								 * therefore, we assume exponential decay of consumption by applying 
								 * the fudgefactor
								 */
								fudgeFactor = 0.0;
								
								for (int  iReac=0;iReac<_reactions.size();iReac++)
								//compute the new reaction
								fudgeFactor += _biomass[iReac]._conc.grid[i][j][k]*
											_reactions.get(iReac).getSoluteYield()[iSolute]*
											singleReac[iReac].grid[i][j][k];
								
								fudgeFactor /= _solute[iSolute]._conc.grid[i][j][k];
								fudgeFactor *= CurrentUptakeStep;
								//rescale the reaction to fit the new uptake
								Double oldReac = _solute[iSolute]._reac.grid[i][j][k];
								_solute[iSolute]._reac.grid[i][j][k] = 
										Math.expm1(fudgeFactor)*_solute[iSolute]._conc.grid[i][j][k]/CurrentUptakeStep;
								
								for (int  iReac=0;iReac<_reactions.size();iReac++)
	
									// Only do this to reactions that consume this solute
									if (_reactions.get(iReac).getSoluteYield()[iSolute] < 0.0)
									{
										singleReac[iReac].grid[i][j][k] *= oldReac/_solute[iSolute]._reac.grid[i][j][k];
										LogFile.writeLogAlways("Reaction '"+_reactions.get(iReac).reactionName+
																"' is getting rescaled");
									}
									
								
									
								if (_soluteIndex.contains(iSolute)) // we do no uptake for steadyState-Solutes
									_solute[iSolute]._conc.grid[i][j][k] *= Math.exp(fudgeFactor);
								
							}
							// the usual behaviour
							else if (_soluteIndex.contains(iSolute)) // we do no uptake for steadyState-Solutes
									_solute[iSolute]._conc.grid[i][j][k]+= soluteDiff; // add/subtract our reacted substrate
						}
			// move our time-counter forward
			time=time+CurrentUptakeStep;
			// finish the uptakeStep by setting our realGrid
			for (int i=0;i<nSolute;i++)
				_solute[i].applyComputation(); //needs to be set, because steadyState always solves from realGrid
		}
	}
	
	/**
	 * runs the implicit solver in the finest order
	 */
	public void stepSolveDiffusionReactionNoMultiGrid() 
	{
		//start the actual solver by setting your concentration (initial or previous timestep)
		for (int iSolute : _soluteIndex)
		{
			_solute[iSolute].setSoluteGridToLastTimeStep(_solute[iSolute]._conc, _bLayer);
		}
		//solves our domain
		coupledSolver();
		
		//write the solved domain to our grid
		for (int iSolute : _soluteIndex)
			_solute[iSolute].applyComputation();
	}	
	
	/**
	 * Update concentration in the reactor
	 */
	public void updateBulk() 
	{
		// Update reaction rates
		// this yields solute change rates in fg.L-1.hr-1
		updateReacRateAndDiffRate();

		// Find the connected bulks and agars and update their concentration
		for (AllBC aBC : myDomain.getAllBoundaries()) 
		{
			if (aBC.hasBulk()) 
				aBC.updateBulk(allSolute, allReac, internTimeStep);
			if (aBC.hasAgar()) 
				aBC.updateAgar(allSolute, allReac, internTimeStep);
		}

		// Refresh the bulk concentration of the multigrids
		for (int iSolute : _soluteIndex)
			_solute[iSolute].readBulk();
	}

	/**
	 * Call all the agents and read their uptake-rate for the current
	 * concentration
	 */
	public void updateReacRateAndDiffRate()
	{
		Double yield, specRate, biomass;
		// Reset rates and derivative rates grids
		for (int iSolute=0; iSolute<nSolute; iSolute++) 
		{ 
			_solute[iSolute].resetReaction();
			allSolute[iSolute] = _solute[iSolute]._conc;
			allReac[iSolute] = _solute[iSolute]._reac;
			allDiffReac[iSolute] = _solute[iSolute]._diffReac; 
		}
		for (int iReac = 0; iReac<_reactions.size(); iReac++)
		{
			singleReac[iReac].setAllValueAt(0.0);
			_reactions.get(iReac).applySingleReaction(allSolute, singleReac[iReac], allDiffReac,
					_biomass[iReac]._conc);
			for ( int iSolute : _soluteIndex )
			{
				yield = _reactions.get(iReac).getSoluteYield()[iSolute];
				if ( yield == 0.0 )
					continue;
				for (int i=0;i<_solute[iSolute]._reac.grid.length;i++)
					for (int j=0;j<_solute[iSolute]._reac.grid[1].length;j++)
						for (int k=0;k<_solute[iSolute]._reac.grid[1][1].length;k++)
						{
							specRate = singleReac[iReac].grid[i][j][k];
							biomass = _biomass[iReac]._conc.grid[i][j][k];
							allReac[iSolute].grid[i][j][k] += specRate * biomass * yield;
						}
			}
		}
	}
}
