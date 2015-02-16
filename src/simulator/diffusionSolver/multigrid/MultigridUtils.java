/**
 * \package diffusionSolver.multigrid
 * \brief Package of classes used to aid solver calculation for multi-grid scenarios.
 * 
 * Package of classes used to capture the diffusion solvers that can be defined in the protocol file. This package is 
 * part of iDynoMiCS v1.2, governed by the CeCILL license under French law and abides by the rules of distribution of free software.  
 * You can use, modify and/ or redistribute iDynoMiCS under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at 
 * the following URL  "http://www.cecill.info".
 */
package simulator.diffusionSolver.multigrid;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import utils.ExtraMath;
import simulator.SoluteGrid;

public abstract class MultigridUtils {

	// boundary layer threshold
	private static final double BLTHRESH  = MultigridSolute.BLTHRESH;
	public static final String  SEPARATOR = " ";

	/**
	 * retuns log2(n - 1) * 0.5
	 * 
	 * @param n
	 * @return order of multigrid
	 * @throws Exception
	 */
	public static int order(Integer n) throws Exception 
	{
		Double out = ExtraMath.log2(n - 1.0);
		if (out%1==0) return out.intValue();
		else throw (new Exception("invalid grid value ("+n+") must be 2^i + 1"));
	}

	/**
	 * \brief Restricts the data in matrix u to a grid one order coarser. 
	 * Restriction excludes border points.
	 * 
	 * Restricts the data in matrix u to a grid one order coarser. 
	 * Restriction excludes border points.
	 * 
	 * TODO currently unused
	 * 
	 * @param fineGrid finer grid
	 * @param coarsegrid coarser grid
	 */
	public static void restrict8c(SoluteGrid fineGrid, SoluteGrid coarsegrid) 
	{
		double[][][] uc = coarsegrid.grid;
		double[][][] u = fineGrid.grid;

		int lc = uc[0][0].length-2;
		int mc = uc[0].length-2;
		int nc = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (lc==1 ? 1.0f/16.0f : 1.0f/32.0f); // pre-compute

		for (k = 1, kc = 1; kc<=lc; kc++, k += 2)
			for (j = 1, jc = 1; jc<=mc; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nc; ic++, i += 2) 
				{
					// 4-connectivity weight
					uc[ic][jc][kc] = 2*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]);
					// 8-connectivity weight
					uc[ic][jc][kc] += u[i+1][j+1][k]+u[i+1][j-1][k]+u[i-1][j+1][k]+u[i-1][j-1][k];
					// 3rd dimension (4-C)
					uc[ic][jc][kc] += 2*(lc==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
					// 3rd dimension (8-C)					
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i-1][j-1][k-1]+u[i-1][j][k-1]+u[i-1][j+1][k-1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j-1][k-1]+u[i][j+1][k-1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i+1][j-1][k-1]+u[i+1][j][k-1]+u[i+1][j+1][k-1]);
					
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i-1][j-1][k+1]+u[i-1][j][k+1]+u[i-1][j+1][k+1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j-1][k+1]+u[i][j+1][k+1]);
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i+1][j-1][k+1]+u[i+1][j][k+1]+u[i+1][j+1][k+1]);
					
					uc[ic][jc][kc] += 4*u[i][j][k];
					uc[ic][jc][kc] *= nfac;
				}
		coarsegrid.refreshBoundary("bLayer");
	}

	/**
	 * same as restrict8c, but simpler
	 * 
	 * @param fineGrid
	 * @param coarseGrid
	 * @param behaviour determines the behaviour of the refresh-boundary, because we use
	 * this restriction for more than just concentration-grids
	 */
	public static void restrict4c(SoluteGrid fineGrid, SoluteGrid coarseGrid, String behaviour) 
	{

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int lc = uc[0][0].length-2;
		int mc = uc[0].length-2;
		int nc = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (lc==1 ? 1.0f/8.0f : 1.0f/12.0f); // pre-compute

		for (k = 1, kc = 1; kc<=lc; kc++, k += 2)
			for (j = 1, jc = 1; jc<=mc; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nc; ic++, i += 2) {
					// special case for 2D (when lc = 1)
					uc[ic][jc][kc] = u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k];
					uc[ic][jc][kc] += (lc==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
					uc[ic][jc][kc] *= nfac;
					// first half of new value is computed
					uc[ic][jc][kc] += 0.5f*u[i][j][k]; // adds the other half
					
					if(Double.isNaN(uc[ic][jc][kc]))
					{
						System.out.print(1);
					}
				}
		coarseGrid.refreshBoundary(behaviour);
	}

	/**
	 * Restricts the data in matrix u to a grid one order coarser. Restriction
	 * excludes border points for points inside the boundary layer, defined by
	 * data in bl. Restriction excludes border points and points outside the
	 * boundary layer (where bl >= 0.5). Points outside boundary layer are
	 * skipped and, therefore, preserve their original value.
	 *  
	 * TODO currently unused
	 * 
	 * @param fineGrid	finer grid
	 * @param coarseGrid	coarser grid
	 * @param bl	boundary layer at coarser grid
	 */
	public static void restrictBoundaryLayer8c(SoluteGrid fineGrid, SoluteGrid coarseGrid,
	        																double[][][] bl) 
	{

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = uc[0][0].length-2;
		int nJ = uc[0].length-2;
		int nI = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		float nfac = (nK==1 ? 1.0f/16.0f : 1.0f/32.0f); // pre-compute

		// full weighting restriction-algorithm
		for (k = 1, kc = 1; kc<=nK; kc++, k += 2)
			for (j = 1, jc = 1; jc<=nJ; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nI; ic++, i += 2)
					if (bl[ic][jc][kc]>=BLTHRESH) 
					{
						// 4-connectivity weight
						uc[ic][jc][kc] = 2*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]);
						// 8-connectivity weight
						uc[ic][jc][kc] += u[i+1][j+1][k]+u[i+1][j-1][k]+u[i-1][j+1][k]+u[i-1][j-1][k];
						// 3rd dimension (4-C)
						uc[ic][jc][kc] += 2*(nK==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
						// 3rd dimension (8-C)					
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i-1][j-1][k-1]+u[i-1][j][k-1]+u[i-1][j+1][k-1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j-1][k-1]+u[i][j+1][k-1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i+1][j-1][k-1]+u[i+1][j][k-1]+u[i+1][j+1][k-1]);
						
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i-1][j-1][k+1]+u[i-1][j][k+1]+u[i-1][j+1][k+1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j-1][k+1]+u[i][j+1][k+1]);
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i+1][j-1][k+1]+u[i+1][j][k+1]+u[i+1][j+1][k+1]);
						
						uc[ic][jc][kc] += 4*u[i][j][k];
						uc[ic][jc][kc] *= nfac;
					}
		coarseGrid.refreshBoundary("conc");
	}

	/**
	 * same as restictboundarylayer8c, but simpler
	 * 
	 * @param fineGrid
	 * @param coarseGrid
	 * @param bl
	 */
	public static void restrictBoundaryLayer4c(SoluteGrid fineGrid, SoluteGrid coarseGrid,
	        																double[][][] bl) 
	{

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = uc[0][0].length-2;
		int nJ = uc[0].length-2;
		int nI = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		Double nfac = (nK==1 ? 1.0/8.0 : 1.0/12.0); // pre-compute

		// half-weighting restriction-algorithm
		for (k = 1, kc = 1; kc<=nK; kc++, k += 2)
			for (j = 1, jc = 1; jc<=nJ; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nI; ic++, i += 2)
					if (bl[ic][jc][kc]>=BLTHRESH) 
					{
						// special case for 2D (when lc = 1)
						uc[ic][jc][kc] = u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k];
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
						uc[ic][jc][kc] *= nfac;
						uc[ic][jc][kc] += 0.5f*u[i][j][k];
					}
		
		coarseGrid.refreshBoundary("conc");
	}	
	
	/**
	 * made for recoarsening in the time-dependent method TODO currently unnecessary
	 * 
	 * @param fineGrid
	 * @param coarseGrid
	 * @param bl	the boundary Layer
	 * @param sBulk
	 */
	public static void restrictStep(SoluteGrid fineGrid, SoluteGrid coarseGrid,
	        										double[][][] bl, double sBulk) 
	{

		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = uc[0][0].length-2;
		int nJ = uc[0].length-2;
		int nI = uc.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// implements 2D and 3D
		Double nfac = (nK==1 ? 1.0/8.0 : 1.0/12.0); // pre-compute

		// half-weighting restriction-algorithm
		for (k = 1, kc = 1; kc<=nK; kc++, k += 2)
			for (j = 1, jc = 1; jc<=nJ; jc++, j += 2)
				for (i = 1, ic = 1; ic<=nI; ic++, i += 2)
					if (bl[ic][jc][kc]>=BLTHRESH) 
					{
						// special case for 2D (when lc = 1)
						uc[ic][jc][kc] = u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k];
						uc[ic][jc][kc] += (nK==1 ? 0.0f : u[i][j][k+1]+u[i][j][k-1]);
						uc[ic][jc][kc] *= nfac;
						uc[ic][jc][kc] += 0.5f*u[i][j][k];
					}
					else
						uc[ic][jc][kc] = sBulk;
		
		coarseGrid.refreshBoundary("conc");
	}
	
	/**
	 * Interpolates the data in matrix uc to a grid one order finer for cubic
	 * matrices. Interpolation excludes border points.
	 * 
	 * @param fineGrid finer grid
	 * @param coarsegrid coarser grid
	 */
	static void interpolate(SoluteGrid fineGrid, SoluteGrid coarsegrid) 
	{
		double[][][] uc = coarsegrid.grid;
		double[][][] u = fineGrid.grid;

		int l = u[0][0].length-2;
		int m = u[0].length-2;
		int n = u.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// copy points
		for (kc = 1, k = 1; k<=l; kc++, k += 2) {
			for (jc = 1, j = 1; j<=m; jc++, j += 2) {
				for (ic = 1, i = 1; i<=n; ic++, i += 2)
					u[i][j][k] = uc[ic][jc][kc];
			}
		}
		// interpolate vertically
		for (k = 1; k<=l; k += 2) {
			for (j = 1; j<=m; j += 2) {
				for (i = 2; i<n; i += 2)
					u[i][j][k] = 0.5f*(u[i+1][j][k]+u[i-1][j][k]);
			}
		}
		// interpolate sideways
		for (k = 1; k<=l; k += 2) {
			for (j = 2; j<m; j += 2) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j+1][k]+u[i][j-1][k]);
			}
		}
		for (k = 2; k<l; k += 2) {
			for (j = 1; j<=m; j++) {
				for (i = 1; i<=n; i++)
					u[i][j][k] = 0.5f*(u[i][j][k+1]+u[i][j][k-1]);
			}
		}

		fineGrid.refreshBoundary("conc");
	}
	
	/**
	 * restriction-method for a cell-centered grid. 
	 * Currently unused, but you might want to add multigrid-behaviour to the solver again
	 * 
	 * @param fineGrid
	 * @param coarseGrid
	 * @author Max Adolph, TU Dresden
	 */
	static void restrictCellCentered(SoluteGrid fineGrid, SoluteGrid coarseGrid) 
	{
		
		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;
		
		int nK = u[0][0].length-2;
		int nJ = u[0].length-2;
		int nI = u.length-2;
		
		int nKc = uc[0][0].length-2;
		int nJc = uc[0].length-2;
		int nIc = uc.length-2;
		
		//int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid
		
		// restrict inner Points
		for (ic=1;ic<nIc+1;ic++)
			for (jc=1;jc<nJc+1;jc++)
				for (kc=1;kc<nKc+1;kc++)
					uc[ic][jc][kc]=(u[(ic-1)*2+1][(jc-1)*2+1][(kc-1)*2+1]+
									u[(ic-1)*2+2][(jc-1)*2+1][(kc-1)*2+1]+
									u[(ic-1)*2+1][(jc-1)*2+2][(kc-1)*2+1]+
									u[(ic-1)*2+2][(jc-1)*2+2][(kc-1)*2+1]+
									u[(ic-1)*2+1][(jc-1)*2+1][(kc-1)*2+2]+
									u[(ic-1)*2+2][(jc-1)*2+1][(kc-1)*2+2]+
									u[(ic-1)*2+1][(jc-1)*2+2][(kc-1)*2+2]+
									u[(ic-1)*2+2][(jc-1)*2+2][(kc-1)*2+2])/8;
				
		//restrict faces
		
		// Z-Faces
		kc = 0;
		for (ic=1;ic<nIc+1;ic++)
			for (jc=1;jc<nJc+1;jc++)
				uc[ic][jc][kc]=(u[(ic-1)*2+1][(jc-1)*2+1][0]+
								u[(ic-1)*2+2][(jc-1)*2+1][0]+
								u[(ic-1)*2+1][(jc-1)*2+2][0]+
								u[(ic-1)*2+2][(jc-1)*2+2][0])/4;
			
		kc = nKc+1;
		for (ic=1;ic<nIc+1;ic++)
			for (jc=1;jc<nJc+1;jc++)
				uc[ic][jc][kc]=(u[(ic-1)*2+1][(jc-1)*2+1][nK+1]+
								u[(ic-1)*2+2][(jc-1)*2+1][nK+1]+
								u[(ic-1)*2+1][(jc-1)*2+2][nK+1]+
								u[(ic-1)*2+2][(jc-1)*2+2][nK+1])/4;
		
		// Y-Faces
		jc = 0;
		for (ic=1;ic<nIc+1;ic++)
			for (kc=1;kc<nKc+1;kc++)
				uc[ic][jc][kc]=(u[(ic-1)*2+1][0][(kc-1)*2+1]+
								u[(ic-1)*2+2][0][(kc-1)*2+1]+
								u[(ic-1)*2+1][0][(kc-1)*2+2]+
								u[(ic-1)*2+2][0][(kc-1)*2+2])/4;

		jc = nJc+1;
		for (ic=1;ic<nIc+1;ic++)
			for (kc=1;kc<nKc+1;kc++)
				uc[ic][jc][kc]=(u[(ic-1)*2+1][nJ+1][(kc-1)*2+1]+
								u[(ic-1)*2+2][nJ+1][(kc-1)*2+1]+
								u[(ic-1)*2+1][nJ+1][(kc-1)*2+2]+
								u[(ic-1)*2+2][nJ+1][(kc-1)*2+2])/4;
		
		// X-Faces
		ic = 0;
		for (jc=1;jc<nJc+1;jc++)
			for (kc=1;kc<nKc+1;kc++)
				ic = 0;
				uc[ic][jc][kc]=(u[0][(jc-1)*2+1][(kc-1)*2+1]+
								u[0][(jc-1)*2+2][(kc-1)*2+1]+
								u[0][(jc-1)*2+1][(kc-1)*2+2]+
								u[0][(jc-1)*2+2][(kc-1)*2+2])/4;
		
		ic = nIc+1;
		for (jc=1;jc<nJc+1;jc++)
			for (kc=1;kc<nKc+1;kc++)
				ic = nIc+1;
				uc[ic][jc][kc]=(u[nI+1][(jc-1)*2+1][(kc-1)*2+1]+
								u[nI+1][(jc-1)*2+2][(kc-1)*2+1]+
								u[nI+1][(jc-1)*2+1][(kc-1)*2+2]+
								u[nI+1][(jc-1)*2+2][(kc-1)*2+2])/4;
			
	}
	
	/**
	 * interpolates the next finer grid in a cell-centered manner. This means, that you can't copy points, 
	 * but instead have to interpolate all of them. just imagine a coarse Grid [__][__][__] and the smaller grid [][][][][][]. 
	 * Then the old points would be 0.5, 1.5, 2.5 and the new points would be at 0.25, 0.75, 1.25, 1.75, 2.25, 2.75 
	 * (assuming 0 is the left border and every element is 1 long) 
	 * 
	 * currently unused, because the time-dependent solver does not use multigrid and the steady-state-Solver still
	 * uses a node-centered grid
	 * 
	 * @param fineGrid
	 * @param coarseGrid
	 * @author Max Adolph, TU Dresden
	 */
	static void interpolateCellCentered(SoluteGrid fineGrid, SoluteGrid coarseGrid) 
	{
		
		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;
		
		int nK = u[0][0].length-2;
		int nJ = u[0].length-2;
		int nI = u.length-2;
		
		int nKc = uc[0][0].length-2;
		int nJc = uc[0].length-2;
		int nIc = uc.length-2;
		
		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid
		
		//3D-case
		if (nK>1)
		{
			for (jc=1;jc<nJc;jc++)
				for (kc=1;kc<nKc;kc++)
					for (ic=1;ic<nIc;ic++)
						// 0 shifted indices = *27; 1 shifted index = *9; 2 shifted indices = *3; 3 shifted indices = *1
						// [the obscure coefficients are the result of some matrix-operations to spare some steps]
						// restored original behaviour => TODO test for correctness
						// imagine a cube of fine points inside the cube of coarse points with 
						// its face towards (face 1) and away (face 2) from you 
					{
						//top    left  of face 1
						u[(ic-1)*2+2][(jc-1)*2+2][(kc-1)*2+2]=(27*uc[ic][jc][kc]+9*(uc[ic+1][jc][kc]+uc[ic][jc+1][kc])+3*uc[ic+1][jc+1][kc]
								+9*uc[ic][jc][kc+1]+3*(uc[ic+1][jc][kc+1]+uc[ic][jc+1][kc+1])+uc[ic+1][jc+1][kc+1])/64;
						
						//top    right of face 1
						u[(ic-1)*2+3][(jc-1)*2+2][(kc-1)*2+2]=(27*uc[ic+1][jc][kc]+9*(uc[ic+1][jc+1][kc]+uc[ic][jc+1][kc])+3*uc[ic][jc+1][kc]
								+9*uc[ic+1][jc][kc+1]+3*(uc[ic][jc][kc+1]+uc[ic+1][jc+1][kc+1])+uc[ic][jc+1][kc+1])/64;
						
						//bottom right of face 1
						u[(ic-1)*2+3][(jc-1)*2+3][(kc-1)*2+2]=(27*uc[ic+1][jc+1][kc]+9*(uc[ic+1][jc][kc]+uc[ic][jc+1][kc])+3*uc[ic][jc][kc]
								+9*uc[ic+1][jc+1][kc+1]+3*(uc[ic+1][jc][kc+1]+uc[ic][jc+1][kc+1])+uc[ic][jc][kc+1])/64;
						
						//bottom left  of face 1
						u[(ic-1)*2+2][(jc-1)*2+3][(kc-1)*2+2]=(27*uc[ic][jc+1][kc]+9*(uc[ic+1][jc+1][kc]+uc[ic][jc][kc])+3*uc[ic+1][jc][kc]
								+9*uc[ic][jc+1][kc+1]+3*(uc[ic+1][jc+1][kc+1]+uc[ic][jc][kc+1])+uc[ic+1][jc][kc+1])/64;
						
						//top    left  of face 2
						u[(ic-1)*2+2][(jc-1)*2+2][(kc-1)*2+3]=(9*uc[ic][jc][kc]+3*(uc[ic+1][jc][kc]+uc[ic][jc+1][kc])+uc[ic+1][jc+1][kc]
								+27*uc[ic][jc][kc+1]+9*(uc[ic+1][jc][kc+1]+uc[ic][jc+1][kc+1])+3*uc[ic+1][jc+1][kc+1])/64;
						
						//top    right of face 2
						u[(ic-1)*2+3][(jc-1)*2+2][(kc-1)*2+3]=(9*uc[ic+1][jc][kc]+3*(uc[ic][jc][kc]+uc[ic+1][jc+1][kc])+uc[ic][jc+1][kc]
								+27*uc[ic+1][jc][kc+1]+9*(uc[ic][jc][kc+1]+uc[ic+1][jc+1][kc+1])+3*uc[ic][jc+1][kc+1])/64;
						
						//bottom right of face 2
						u[(ic-1)*2+3][(jc-1)*2+3][(kc-1)*2+3]=(9*uc[ic+1][jc+1][kc]+3*(uc[ic+1][jc][kc]+uc[ic][jc+1][kc])+uc[ic][jc][kc]
								+27*uc[ic+1][jc+1][kc+1]+9*(uc[ic+1][jc][kc+1]+uc[ic][jc+1][kc+1])+3*uc[ic][jc][kc+1])/64;
						
						//bottom left  of face 2
						u[(ic-1)*2+2][(jc-1)*2+3][(kc-1)*2+3]=(9*uc[ic][jc+1][kc]+3*(uc[ic+1][jc+1][kc]+uc[ic][jc][kc])+uc[ic+1][jc][kc]
								+27*uc[ic][jc+1][kc+1]+9*(uc[ic+1][jc+1][kc+1]+uc[ic][jc][kc+1])+3*uc[ic+1][jc][kc+1])/64;
					}	

			/**
			 * that was the actual work, but we need to cope with the fact that our domain has
			 * no corners and edges... so here comes the ugly part
			 */
			
			//set edges (just to be complete)
			for (k=0;k<nK+2;k++)
				for (j=0;j<nJ+2;j++)
					for (i=0;i<nI+2;i=i+nI+1)
						u[i][j][k]=0;

			for (j=0;j<nJ+2;j++)
				for (i=0;i<nI+2;i++)
					for (k=0;k<nK+2;k=k+nK+1)
						u[i][j][k]=0;

			for (i=0;i<nI+2;i++)
				for (k=0;k<nK+2;k++)
					for (j=0;j<nJ+2;j=j+nJ+1)
						u[i][j][k]=0;

			//initiate fixed values
			int m,n;

			//set faces of padding
			m=nJc+1;
			for (ic=1;ic<nIc;ic++)
				for (kc=1;kc<nKc;kc++)
				{
					u[(ic-1)*2+2][0][(kc-1)*2+2]=(9*uc[ic][0][kc]+3*(uc[ic+1][0][kc]+uc[ic][0][kc+1])+uc[ic+1][0][kc+1])/16;
					u[(ic-1)*2+3][0][(kc-1)*2+2]=(3*(uc[ic][0][kc]+uc[ic+1][0][kc+1])+9*uc[ic+1][0][kc]+uc[ic][0][kc+1])/16;
					u[(ic-1)*2+2][0][(kc-1)*2+3]=(3*(uc[ic+1][0][kc+1]+uc[ic][0][kc])+9*uc[ic][0][kc+1]+uc[ic+1][0][kc])/16;
					u[(ic-1)*2+3][0][(kc-1)*2+3]=(9*uc[ic+1][0][kc+1]+3*(uc[ic][0][kc+1]+uc[ic+1][0][kc])+uc[ic][0][kc])/16;

					u[(ic-1)*2+2][nJ+1][(kc-1)*2+2]=(9*uc[ic][m][kc]+3*(uc[ic+1][m][kc]+uc[ic][m][kc+1])+uc[ic+1][m][kc+1])/16;
					u[(ic-1)*2+3][nJ+1][(kc-1)*2+2]=(3*(uc[ic][m][kc]+uc[ic+1][m][kc+1])+9*uc[ic+1][m][kc]+uc[ic][m][kc+1])/16;
					u[(ic-1)*2+2][nJ+1][(kc-1)*2+3]=(3*(uc[ic+1][m][kc+1]+uc[ic][m][kc])+9*uc[ic][m][kc+1]+uc[ic+1][m][kc])/16;
					u[(ic-1)*2+3][nJ+1][(kc-1)*2+3]=(9*uc[ic+1][m][kc+1]+3*(uc[ic][m][kc+1]+uc[ic+1][m][kc])+uc[ic][m][kc])/16;

				}
			m=nIc+1;
			for (jc=1;jc<nJc;jc++)
				for (kc=1;kc<nKc;kc++)
				{

					u[0][(jc-1)*2+2][(kc-1)*2+2]=(9*uc[0][jc][kc]+3*(uc[0][jc+1][kc]+uc[0][jc][kc+1])+uc[0][jc+1][kc+1])/16;
					u[0][(jc-1)*2+3][(kc-1)*2+2]=(3*(uc[0][jc][kc]+uc[0][jc+1][kc+1])+9*uc[0][jc+1][kc]+uc[0][jc][kc+1])/16;
					u[0][(jc-1)*2+2][(kc-1)*2+3]=(3*(uc[0][jc+1][kc+1]+uc[0][jc][kc])+9*uc[0][jc][kc+1]+uc[0][jc+1][kc])/16;
					u[0][(jc-1)*2+3][(kc-1)*2+3]=(9*uc[0][jc+1][kc+1]+3*(uc[0][jc][kc+1]+uc[0][jc+1][kc])+uc[0][jc][kc])/16;

					u[nI+1][(jc-1)*2+2][(kc-1)*2+2]=(9*uc[m][jc][kc]+3*(uc[m][jc+1][kc]+uc[m][jc][kc+1])+uc[m][jc+1][kc+1])/16;
					u[nI+1][(jc-1)*2+3][(kc-1)*2+2]=(3*(uc[m][jc][kc]+uc[m][jc+1][kc+1])+9*uc[m][jc+1][kc]+uc[m][jc][kc+1])/16;
					u[nI+1][(jc-1)*2+2][(kc-1)*2+3]=(3*(uc[m][jc+1][kc+1]+uc[m][jc][kc])+9*uc[m][jc][kc+1]+uc[m][jc+1][kc])/16;
					u[nI+1][(jc-1)*2+3][(kc-1)*2+3]=(9*uc[m][jc+1][kc+1]+3*(uc[m][jc][kc+1]+uc[m][jc+1][kc])+uc[m][jc][kc])/16;

				}
			m=nKc+1;
			for (ic=1;ic<nIc;ic++)
				for (jc=1;jc<nJc;jc++)
				{

					u[(ic-1)*2+2][(jc-1)*2+2][0]=(9*uc[ic][jc][0]+3*(uc[ic+1][jc][0]+uc[ic][jc+1][0])+uc[ic+1][jc+1][0])/16;
					u[(ic-1)*2+3][(jc-1)*2+2][0]=(3*(uc[ic][jc][0]+uc[ic+1][jc+1][0])+9*uc[ic+1][jc][0]+uc[ic][jc+1][0])/16;
					u[(ic-1)*2+2][(jc-1)*2+3][0]=(3*(uc[ic+1][jc+1][0]+uc[ic][jc][0])+9*uc[ic][jc+1][0]+uc[ic+1][jc][0])/16;
					u[(ic-1)*2+3][(jc-1)*2+3][0]=(9*uc[ic+1][jc+1][0]+3*(uc[ic][jc+1][0]+uc[ic+1][jc][0])+uc[ic][jc][0])/16;

					u[(ic-1)*2+2][(jc-1)*2+2][nK+1]=(9*uc[ic][jc][m]+3*(uc[ic+1][jc][m]+uc[ic][jc+1][m])+uc[ic+1][jc+1][m])/16;
					u[(ic-1)*2+3][(jc-1)*2+2][nK+1]=(3*(uc[ic][jc][m]+uc[ic+1][jc+1][m])+9*uc[ic+1][jc][m]+uc[ic][jc+1][m])/16;
					u[(ic-1)*2+2][(jc-1)*2+3][nK+1]=(3*(uc[ic+1][jc+1][m]+uc[ic][jc][m])+9*uc[ic][jc+1][m]+uc[ic+1][jc][m])/16;
					u[(ic-1)*2+3][(jc-1)*2+3][nK+1]=(9*uc[ic+1][jc+1][m]+3*(uc[ic][jc+1][m]+uc[ic+1][jc][m])+uc[ic][jc][m])/16;

				}

			//set edges of padding
			m=nJc+1;
			n=nIc+1;
			for (kc=1;kc<nKc;kc++)
			{
				u[1][0][(kc-1)*2+2]=(3*uc[1][0][kc]+uc[1][0][kc+1])/4.0;
				u[1][0][(kc-1)*2+3]=(3*uc[1][0][kc+1]+uc[1][0][kc])/4.0;
				u[0][1][(kc-1)*2+2]=(3*uc[0][1][kc]+uc[0][1][kc+1])/4.0;
				u[0][1][(kc-1)*2+3]=(3*uc[0][1][kc+1]+uc[0][1][kc])/4.0;

				u[1][nJ+1][(kc-1)*2+2]=(3*uc[1][m][kc]+uc[1][m][kc+1])/4.0;
				u[1][nJ+1][(kc-1)*2+3]=(3*uc[1][m][kc+1]+uc[1][m][kc])/4.0;
				u[0][nJ][(kc-1)*2+2]=(3*uc[0][m-1][kc]+uc[0][m-1][kc+1])/4.0;
				u[0][nJ][(kc-1)*2+3]=(3*uc[0][m-1][kc+1]+uc[0][m-1][kc])/4.0;

				u[nI+1][1][(kc-1)*2+2]=(3*uc[n][0][kc]+uc[n][1][kc+1])/4.0;
				u[nI+1][1][(kc-1)*2+3]=(3*uc[n][0][kc+1]+uc[n][1][kc])/4.0;
				u[nI][0][(kc-1)*2+2]=(3*uc[n-1][0][kc]+uc[n-1][0][kc+1])/4.0;
				u[nI][0][(kc-1)*2+3]=(3*uc[n-1][0][kc+1]+uc[n-1][0][kc])/4.0;

				u[nI][nJ+1][(kc-1)*2+2]=(3*uc[n-1][m][kc]+uc[n-1][m][kc+1])/4.0;
				u[nI][nJ+1][(kc-1)*2+3]=(3*uc[n-1][m][kc+1]+uc[n-1][m][kc])/4.0;
				u[nI+1][nJ][(kc-1)*2+2]=(3*uc[n][m-1][kc]+uc[n][m-1][kc+1])/4.0;
				u[nI+1][nJ][(kc-1)*2+3]=(3*uc[n][m-1][kc+1]+uc[n][m-1][kc])/4.0;
			}
			m=nIc+1;
			n=nKc+1;
			for (jc=1;jc<nJc;jc++)
			{
				u[1][(jc-1)*2+2][0]=(3*uc[1][jc][0]+uc[1][jc+1][0])/4.0;
				u[1][(jc-1)*2+3][0]=(3*uc[1][jc+1][0]+uc[1][jc][0])/4.0;
				u[0][(jc-1)*2+2][1]=(3*uc[0][jc][1]+uc[0][jc+1][1])/4.0;
				u[0][(jc-1)*2+3][1]=(3*uc[0][jc+1][1]+uc[0][jc][1])/4.0;

				u[nI+1][(jc-1)*2+2][1]=(3*uc[m][jc][1]+uc[m][jc+1][1])/4.0;
				u[nI+1][(jc-1)*2+3][1]=(3*uc[m][jc+1][1]+uc[m][jc][1])/4.0;
				u[nI][(jc-1)*2+2][0]=(3*uc[m-1][jc][0]+uc[m-1][jc+1][0])/4.0;
				u[nI][(jc-1)*2+3][0]=(3*uc[m-1][jc+1][0]+uc[m-1][jc][0])/4.0;

				u[1][(jc-1)*2+2][nK+1]=(3*uc[1][jc][n]+uc[1][jc+1][n])/4.0;
				u[1][(jc-1)*2+3][nK+1]=(3*uc[1][jc+1][n]+uc[1][jc][n])/4.0;
				u[0][(jc-1)*2+2][nK]=(3*uc[0][jc][n-1]+uc[0][jc+1][n-1])/4.0;
				u[0][(jc-1)*2+3][nK]=(3*uc[0][jc+1][n-1]+uc[0][jc][n-1])/4.0;

				u[nI][(jc-1)*2+2][nK+1]=(3*uc[m-1][jc][n]+uc[m-1][jc+1][n])/4.0;
				u[nI][(jc-1)*2+3][nK+1]=(3*uc[m-1][jc+1][n]+uc[m-1][jc][n])/4.0;
				u[nI+1][(jc-1)*2+2][nK]=(3*uc[m][jc][n-1]+uc[m][jc+1][n-1])/4.0;
				u[nI+1][(jc-1)*2+3][nK]=(3*uc[m][jc+1][n-1]+uc[m][jc][n-1])/4.0;
			}
			m=nJc+1;
			n=nKc+1;
			for (ic=1;ic<nIc;ic++)
			{
				u[(ic-1)*2+2][1][0]=(3*uc[ic][1][0]+uc[ic+1][1][0])/4.0;
				u[(ic-1)*2+3][1][0]=(3*uc[ic+1][1][0]+uc[ic][1][0])/4.0;
				u[(ic-1)*2+2][0][1]=(3*uc[ic][0][1]+uc[ic+1][0][1])/4.0;
				u[(ic-1)*2+3][0][1]=(3*uc[ic+1][0][1]+uc[ic][0][1])/4.0;

				u[(ic-1)*2+2][nJ+1][1]=(3*uc[ic][m][1]+uc[ic+1][m][1])/4.0;
				u[(ic-1)*2+3][nJ+1][1]=(3*uc[ic+1][m][1]+uc[ic][m][1])/4.0;
				u[(ic-1)*2+2][nJ][0]=(3*uc[ic][m-1][0]+uc[ic+1][m-1][0])/4.0;
				u[(ic-1)*2+3][nJ][0]=(3*uc[ic+1][m-1][0]+uc[ic][m-1][0])/4.0;

				u[(ic-1)*2+2][1][nK+1]=(3*uc[ic][1][n]+uc[ic+1][1][n])/4.0;
				u[(ic-1)*2+3][1][nK+1]=(3*uc[ic+1][1][n]+uc[ic][1][n])/4.0;
				u[(ic-1)*2+2][0][nK]=(3*uc[ic][0][n-1]+uc[ic+1][0][n-1])/4.0;
				u[(ic-1)*2+3][0][nK]=(3*uc[ic+1][0][n-1]+uc[ic][0][n-1])/4.0;

				u[(ic-1)*2+2][nJ][nK+1]=(3*uc[ic][m-1][n]+uc[ic+1][m-1][n])/4.0;
				u[(ic-1)*2+3][nJ][nK+1]=(3*uc[ic+1][m-1][n]+uc[ic][m-1][n])/4.0;
				u[(ic-1)*2+2][nJ+1][nK]=(3*uc[ic][m][n-1]+uc[ic+1][m][n-1])/4.0;
				u[(ic-1)*2+3][nJ+1][nK]=(3*uc[ic+1][m][n-1]+uc[ic][m][n-1])/4.0;
			}

			//set corners of padding
			for (i=1;i<nI+1;i=i+nI-1)
				for (j=1;j<nJ+1;j=j+nJ-1)
					for (k=0;k<nK+2;k=k+nK+1)
					{
						if (i==nI)
							ic=nIc;
						else
							ic=1;
						if (j==nJ)
							jc=nJc;
						else
							jc=1;
						if (k==nK+1)
							kc=nKc+1;
						else
							kc=0;
						u[i][j][k]=uc[ic][jc][kc];
					}

			for (i=1;i<nI+1;i=i+nI-1)
				for (k=1;k<nK+1;k=k+nK-1)
					for (j=0;j<nJ+2;j=j+nJ+1)
					{
						if (i==nI)
							ic=nIc;
						else
							ic=1;
						if (k==nK)
							kc=nKc;
						else
							kc=1;
						if (j==nJ+1)
							jc=nJc+1;
						else
							jc=0;
						u[i][j][k]=uc[ic][jc][kc];
					}

			for (k=1;k<nK+1;k=k+nK-1)
				for (j=1;j<nJ+1;j=j+nJ-1)
					for (i=0;i<nI+2;i=i+nI+1)
					{
						if (k==nK)
							kc=nKc;
						else
							kc=1;
						if (j==nJ)
							jc=nJc;
						else
							jc=1;
						if (i==nI+1)
							ic=nIc+1;
						else
							ic=0;
						u[i][j][k]=uc[ic][jc][kc];
					}

			//interpolate inner faces
			for (i=2;i<nI;i++)
				for (k=2;k<nK;k++)
				{
					u[i][1][k]=u[i][0][k];
					u[i][nJ][k]=u[i][nJ+1][k];
				}
			for (i=2;i<nI;i++)
				for (j=2;j<nJ;j++)
				{
					u[i][j][1]=u[i][j][0];
					u[i][j][nK]=u[i][j][nK+1];
				}
			for (k=2;k<nK;k++)
				for (j=2;j<nJ;j++)
				{
					u[1][j][k]=u[0][j][k];
					u[nI][j][k]=u[nI+1][j][k];
				}

			for (k=2;k<nK;k++)
			{
				u[1][1][k]=(u[0][1][k]+u[1][0][k])/2;
				u[1][nJ][k]=(u[0][nJ][k]+u[1][nJ+1][k])/2;
				u[nI][1][k]=(u[nI+1][1][k]+u[nI][0][k])/2;
				u[nI][nJ][k]=(u[nI+1][nJ][k]+u[nI][nJ+1][k])/2;
			}
			
			for (j=2;j<nJ;j++)
			{
				u[1][j][1]=(u[0][j][1]+u[1][j][0])/2;
				u[1][j][nK]=(u[0][j][nK]+u[1][j][nK+1])/2;
				u[nI][j][1]=(u[nI+1][j][1]+u[nI][j][0])/2;
				u[nI][nJ][k]=(u[nI+1][j][nK]+u[nI][j][nK+1])/2;
			}
			
			for (i=2;i<nI;i++)
			{
				u[i][1][1]=(u[i][1][0]+u[i][0][1])/2;
				u[i][nJ][1]=(u[i][nJ+1][1]+u[i][nJ][0])/2;
				u[i][1][nK]=(u[i][0][nK]+u[i][1][nK+1])/2;
				u[i][nJ][nK]=(u[i][nJ+1][nK]+u[i][nJ][nK+1])/2;
			}
			
			// set inner corners
			u[1][1][1]=		(u[0][1][1]+		u[1][0][1]+			u[1][1][0])/3;
			u[nI][1][1]=	(u[nI+1][1][1]+		u[nI][0][1]+		u[nI][1][0])/3;
			u[1][nJ][1]=	(u[0][nJ][1]+		u[1][nJ+1][1]+		u[1][nJ][0])/3;
			u[nI][nJ][1]=	(u[nI+1][nJ][1]+	u[nI][nJ+1][1]+		u[nI][nJ][0])/3;
			u[1][1][nK]=	(u[0][1][nK]+		u[1][0][nK]+		u[1][1][nK+1])/3;
			u[nI][1][nK]=	(u[nI+1][1][nK]+	u[nI][0][nK]+		u[nI][1][nK+1])/3;
			u[1][nJ][nK]=	(u[0][nJ][nK]+		u[1][nJ+1][nK]+		u[1][nJ][nK+1])/3;
			u[nI][nJ][nK]=	(u[nI+1][nJ][nK]+	u[nI][nJ+1][nK]+	u[nI][nJ][nK+1])/3;
		}
		
		//2D-case
		if (nK==1)
			for (jc=1;jc<nJc-1;jc++)
				for (ic=1;ic<nIc-1;ic++)
				{
					u[(ic-1)*2+1][(jc-1)*2+1][1]=(9*uc[ic][jc][1]+3*(uc[ic+1][jc][1]+uc[ic][jc+1][1])+uc[ic+1][jc+1][1])/16;
					u[(ic-1)*2+2][(jc-1)*2+1][1]=(3*(uc[ic][jc][1]+uc[ic+1][jc+1][1])+9*uc[ic+1][jc][1]+uc[ic][jc+1][1])/16;
					u[(ic-1)*2+2][(jc-1)*2+2][1]=(3*(uc[ic+1][jc][1]+uc[ic][jc+1][1])+9*uc[ic+1][jc+1][1]+uc[ic][jc][1])/16;
					u[(ic-1)*2+1][(jc-1)*2+2][1]=(9*uc[ic][jc+1][1]+3*(uc[ic][jc][1]+uc[ic+1][jc+1][1])+uc[ic+1][jc][1])/16;
				}
	}
	
	/**
	 * Interpolates the data in matrix uc to a grid one order finer for cubic
	 * matrices for points inside the boundary layer, defined by data in bl.
	 * Interpolation excludes border points and points outside the boundary
	 * layer (where bl >= 0.5). Points outside boundary layer are skipped and,
	 * therefore, preserve their original value.
	 * 
	 * @param fineGrid finer grid
	 * @param coarseGrid coarser grid
	 * @param bl boundary layer at finer grid
	 */
	static void interpolateBoundaryLayer(SoluteGrid fineGrid, SoluteGrid coarseGrid, double[][][] bl) 
	{
		double[][][] uc = coarseGrid.grid;
		double[][][] u = fineGrid.grid;

		int nK = u[0][0].length-2;
		int nJ = u[0].length-2;
		int nI = u.length-2;

		int i, j, k; // indexes for fine grid
		int ic, jc, kc; // indexes for coarse grid

		// copy points
		for (kc = 1, k = 1; k<=nK; kc++, k += 2) 
			for (jc = 1, j = 1; j<=nJ; jc++, j += 2) 
				for (ic = 1, i = 1; i<=nI; ic++, i += 2)
					if (bl[i][j][k]>=BLTHRESH) 
						u[i][j][k] = uc[ic][jc][kc];
			
		// interpolate verically
		for (k = 1; k<=nK; k += 2) 
			for (j = 1; j<=nJ; j += 2) 
				for (i = 2; i<nI; i += 2)
					if (bl[i][j][k]>=BLTHRESH) 
						u[i][j][k] = 0.5*(u[i+1][j][k]+u[i-1][j][k]);
			
		// interpolate sideways
		for (k = 1; k<=nK; k += 2) 
			for (j = 2; j<nJ; j += 2) 
				for (i = 1; i<=nI; i++)
					if (bl[i][j][k]>=BLTHRESH) 
						u[i][j][k] = 0.5*(u[i][j+1][k]+u[i][j-1][k]);
			
		// interpolate in depth
		for (k = 2; k<nK; k += 2) 
			for (j = 1; j<=nJ; j++) 
				for (i = 1; i<=nI; i++)
					if (bl[i][j][k]>=BLTHRESH) 
						u[i][j][k] = 0.5*(u[i][j][k+1]+u[i][j][k-1]);
			
		fineGrid.refreshBoundary("conc");
		
	}

	/**
	 * Set all entries of a matrix to value val
	 * 
	 * @param u
	 * @param val
	 */
	public static void setValues(double u[][][], double val) 
	{
		for (int i = 0; i<u.length; i++)
			for (int j = 0; j<u[i].length; j++)
				for (int k = 0; k<u[i][j].length; k++)
					u[i][j][k] = val;
	}

	/**
	 * Set all entries of a boolean matrix to value val
	 * 
	 * @param u
	 * @param val
	 */
	public static void setValues(boolean u[][][], boolean val)
	{
		for (int i = 0; i<u.length; i++)
			for (int j = 0; j<u[i].length; j++)
				for (int k = 0; k<u[i][j].length; k++)
					u[i][j][k] = val;
	}

	/**
	 * Add every entry of matrix b to the corresponding entry in matrix a
	 * 
	 * @param a
	 * @param b
	 */
	static void addTo(double a[][][], double b[][][]) 
	{
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] += b[i][j][k];
	}

	/**
	 * Subtract every entry of matrix b to the corresponding entry in matrix a.
	 * 
	 * @param a
	 * @param b
	 */
	static void subtractTo(double a[][][], double b[][][]) 
	{
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					a[i][j][k] -= b[i][j][k];
	}

	/**
	 * Create matrix c = a - b
	 * 
	 * @param a
	 * @param b
	 * @return c = a-b
	 */
	public static float[][][] subtract(float a[][][], float b[][][]) 
	{
		int l = a.length;
		int m = a[0].length;
		int n = a[0][0].length;
		float[][][] c = new float[l][m][n];
		for (int i = 0; i<l; i++)
			for (int j = 0; j<m; j++)
				for (int k = 0; k<n; k++)
					c[i][j][k] = a[i][j][k]-b[i][j][k];
		return c;
	}

	/**
	 * Find minimum value in a 3D matrix
	 * 
	 * @param a
	 * @return the minimum value in the matrix
	 */
	public static float min(float a[][][]) 
	{
		float min = a[0][0][0];
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					min = (a[i][j][k]<min ? a[i][j][k] : min);
		return min;
	}

	/**
	 * Find maximum value in a 3D matrix
	 * 
	 * @param a
	 * @return the maximum value in the matrix
	 */
	public static float max(float a[][][]) 
	{
		float max = a[0][0][0];
		for (int i = 0; i<a.length; i++)
			for (int j = 0; j<a[i].length; j++)
				for (int k = 0; k<a[i][j].length; k++)
					max = (a[i][j][k]>max ? a[i][j][k] : max);
		return max;
	}

	/**
	 * compute the euclidian norm of matrix (exceptuating padding)
	 * 
	 * @param a
	 * @return the norm of the matrix
	 */
	public static Double computeNorm(double[][][] a) 
	{
		Double norm = 0.0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					norm += ExtraMath.sq(a[i][j][k]);
		return Math.sqrt(norm);
	}

	/**
	 * @param a
	 * @return the sum of all elements of a
	 */
	public static float computeSum(float[][][] a) 
	{
		float sum = 0;
		for (int i = 1; i<a.length-1; i++)
			for (int j = 1; j<a[i].length-1; j++)
				for (int k = 1; k<a[i][j].length-1; k++)
					sum += a[i][j][k];
		return sum;
	}

	/**
	 * Return values in a matrix (excluding boundaries) as a formatted string
	 * 
	 * @param matrix to output as string
	 * @return string output
	 */
	public static String coreMatrixToString(float[][][] matrix) 
	{
		int n = matrix.length-2;
		int m = matrix[0].length-2;
		int l = matrix[0][0].length-2;
		StringBuffer out = new StringBuffer();
		for (int k = 1; k<=l; k++) {
			for (int i = n; i>=1; i--) {
				for (int j = 1; j<=m; j++) {
					out.append(matrix[i][j][k]);
					// change here for format (presently space separated values
					out.append(SEPARATOR);
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	/**
	 * Return values in a matrix (excluding boundaries) as a formatted string.
	 * This method is used for boolean matrices. Values in output are 1 (for
	 * true) or 0 (for false)
	 * 
	 * @param matrix to output as string
	 * @return string output
	 */
	public static String coreMatrixToString(boolean[][][] matrix) {
		int n = matrix.length-2;
		int m = matrix[0].length-2;
		int l = matrix[0][0].length-2;
		StringBuffer out = new StringBuffer();
		for (int k = 1; k<=l; k++) {
			for (int i = n; i>=1; i--) {
				for (int j = 1; j<=m; j++) {
					out.append(matrix[i][j][k] ? 1 : 0);
					// change here for format (presently space separated values
					out.append(SEPARATOR);
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	/**
	 * Write the full matrix to a string
	 * 
	 * @param matrix
	 * @return a string with the matrix (space separated values)
	 */
	public static String matrixToString(float[][][] matrix) {
		StringBuffer out = new StringBuffer();
		for (int k = 0; k<matrix[0][0].length; k++) {
			for (int i = matrix.length-1; i>=0; i--) {
				for (int j = 0; j<matrix[0].length; j++) {
					out.append(matrix[i][j][k]);
					// change here for format (presently space separated values
					out.append(SEPARATOR);
				}
				out.append("\n");
			}
			out.append("\n");
		}
		return out.toString();
	}

	/**
	 * Create a 2D graphics
	 * 
	 * @param fileName the file to parse
	 * @return 2D matrix
	 */
	public static float[][] readSquareMatrixFromFile(String fileName) {
		String line;
		String[] tokens;
		float[][] dataRead = { { 0 } }; // initialize
		int lineCount = 0;
		int i = 0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while ((line = br.readLine())!=null) {
				br.close();
				tokens = line.split(SEPARATOR);
				int l = tokens.length;
				// in first iteration initialize the 2D matrix
				if (dataRead.length==1) {
					// create the matrix and
					lineCount = l;
					dataRead = new float[lineCount][lineCount];
					i = lineCount-1;
				} else if ((lineCount!=l)&(i>=0)) {

				} else if (i<-1) {

				}
				if (i>=0) {
					// parse the data in the line into a matrix
					for (int j = 0; j<lineCount; j++) {
						dataRead[i][j] = Float.parseFloat(tokens[j]);
					}
				}
				// decrement the line number counter
				i--;
			}
		} catch (FileNotFoundException e) {

		} catch (IOException e) {

		}
		return dataRead;
	}

}