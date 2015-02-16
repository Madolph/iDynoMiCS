/**
 * Project iDynoMiCS (copyright -> see Idynomics.java)
 * ______________________________________________________
 * Class for containing chemical solutes, that are represented by a grid. The
 * grid is padded, 3D grid
 * Diffusivity is expressed in the local time unit
 */

/**
 * @since September 2011
 * @version 1.0
 * @author Katrin Bohl, FSU Jena
 */


package simulator;

import utils.XMLParser;
import utils.UnitConverter;
import simulator.geometry.*;
import simulator.geometry.boundaryConditions.AllBC;

public class SoluteGridTimeDependent extends SpatialGrid {
	// Serial version used for the serialisation of the class
	private static final long serialVersionUID = 1L;
	
	/* ______________________ PROPERTIES OF THE SOLUTE ______________________ */

	// Identifier
	public int                  soluteIndex;
	
	// Diffusivity in water
	public double              diffusivity;
	
	// initial concentration
	public double 			initialConc;
	
	
	// Description of diffusion, carrier and bulk domains for this solute
	private Domain _domain;

	/* _____________________ CONSTRUCTOR ___________________________________ */

	/**
     * Constructor based on the XML file
     * @param aSim
     * @param xmlRoot
     */
	public SoluteGridTimeDependent(Simulator aSim, XMLParser xmlRoot) 
	{
		
		StringBuffer unit = new StringBuffer("");

		gridName = xmlRoot.getAttribute("name");
		soluteIndex = aSim.getSoluteIndex(gridName);
		_domain = aSim.world.getDomain(xmlRoot.getAttribute("domain"));

		/* Set the resolution and create the grid ________________ */
		double resolution = xmlRoot.getParamLength("resolution");
		if (Double.isNaN(resolution)) 
			useDomaingrid();
		else 
			specifyResolution(resolution);
		
		initGrids();

		/* Set the diffusivity ____________________________________ */
		double diffusivityUnit = xmlRoot.getParamDbl("diffusivity", unit);
		diffusivityUnit *= UnitConverter.time(unit.toString());
		diffusivityUnit *= UnitConverter.length(unit.toString());
		diffusivityUnit *= UnitConverter.length(unit.toString());
		diffusivity = diffusivityUnit;

		/* Set the initial concentration __________________________ */
		double concentration = xmlRoot.getParamDbl("concentration");
		// Katrin: If no value specified, set 0.0
		if (Double.isNaN(concentration)) 
			concentration = 0.0;
		setAllValueAt(concentration);
		initialConc=concentration;
		System.out.println("concentration: "+concentration);
	}


	/* The 2 next functions are used when creating the multigrids */

	/**
	 * uses the superclass-method to create this SoluteGrid
	 *  
	 * @param nI	1st dimension
	 * @param nJ	2nd dimension
	 * @param nK	3rd dimension
	 * @param res	resolution
	 */
	public SoluteGridTimeDependent(int nI, int nJ, int nK, double res) 
	{
		super(nI, nJ, nK, res);
	}
	
	/**
	 * uses superclass-method, but adds a name and a domain
	 * 
	 * @param nI	1st dimension
	 * @param nJ	2nd dimension
	 * @param nK	3rd dimension
	 * @param res	resolution
	 * @param aName		name
	 * @param aDomain	domain
	 */
	public SoluteGridTimeDependent(int nI, int nJ, int nK, double res,String aName, Domain aDomain) 
	{
		super(nI, nJ, nK, res);
		gridName = aName;
		_domain = aDomain;
	}
	
	/**
	 * uses an external SoluteGrid to initialise the new SoluteGrid
	 * 
	 * @param nI	1st dimension
	 * @param nJ	2nd dimension
	 * @param nK	3rd dimension
	 * @param res	resolution
	 * @param aSolG	the external grid
	 */
	public SoluteGridTimeDependent(int nI, int nJ, int nK, double res,SoluteGridTimeDependent aSolG) 
	{
		super(nI, nJ, nK, res);
		useExternalSoluteGrid(aSolG);
	}
	
	/**
	 * uses just an external SoluteGrid to initialise a new one
	 * 
	 * @param aSolG the soluteGrid to be mirrored
	 */
	public SoluteGridTimeDependent(SoluteGridTimeDependent aSolG)
	{
		gridName = aSolG.gridName;
		diffusivity = aSolG.diffusivity;
		_domain = aSolG._domain;
		
		_reso = aSolG.getResolution();
		_nI = aSolG.getGridSizeI();
		_nJ = aSolG.getGridSizeJ();
		_nK = aSolG.getGridSizeK();
		
		initGrids();
		
	}

	/**
	 * used in the creation of a new SoluteGrid via an existing SoluteGrid
	 * 
	 * @param aSolG	the SoluteGrid that donates its values
	 */
	public void useExternalSoluteGrid(SoluteGridTimeDependent aSolG) {
		gridName = aSolG.gridName;
		soluteIndex = aSolG.soluteIndex;
		diffusivity = aSolG.diffusivity;
		_domain = aSolG._domain;
	}

	
	/**
     * Use the size and the resolution used to define the computation domain to
     * define the solute grid
     */
	public void useDomaingrid() {
		_reso = _domain.getGrid().getResolution();
		_nI = _domain.getGrid().getGridSizeI();
		System.out.println("i: "+_nI);
		_nJ = _domain.getGrid().getGridSizeJ();
		System.out.println("j: "+_nJ);
		_nK = _domain.getGrid().getGridSizeK();
		System.out.println("k: "+_nK);
	}

	/**
     * Give size of grid for the given resolution, based on length defined in
     * the domain
     * 
     * @param reso resolution
     */
	public void specifyResolution(double reso) {
		_reso = reso;
		_nI = (int) Math.ceil(_domain.getGrid().getGridLength(1)/_reso);
		System.out.println("i: "+_nI);
		_nJ = (int) Math.ceil(_domain.getGrid().getGridLength(2)/_reso);
		System.out.println("j: "+_nJ);
		_nK = (int) Math.ceil(_domain.getGrid().getGridLength(3)/_reso);
		System.out.println("k: "+_nK);
	}

	/* ________________________ MAIN METHODS ______________________________ */

	/**
	 * refreshes the padding of the grid according to the boundaries of the domain
	 * 
	 * @param type the behaviour of the boundary (conc, relDiff or bLayer)
	 */
	public void refreshBoundary(String type) 
	{
		
		for (AllBC aBC:_domain.getAllBoundaries()) 
			aBC.refreshBoundary(this, type);
	}

	/* ________________________ GET & SET __________________________________ */

	/**
	 * returns the name of a grid
	 * 
	 * @return the gridName
	 */
	public String getName() 
	{
		return gridName;
	}

	/**
	 * returns the relative diffusivity-grid of the domain
	 * 
	 * @return the relative diffusivity-grid
	 */
	public double getDiffusivity() 
	{
		return diffusivity;
	}

	/**
	 * returns the domain
	 * 
	 * @return the domain
	 */
	public Domain getDomain() 
	{
		return _domain;
	}

}
