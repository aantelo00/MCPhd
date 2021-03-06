/*
 *  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 2004-2007  Miguel Rojas <miguel.rojas@uni-koeln.de>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.qsar.descriptors.bond;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.charges.GasteigerPEPEPartialCharges;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.qsar.AbstractBondDescriptor;
import org.openscience.cdk.qsar.DescriptorSpecification;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.tools.LonePairElectronChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.Iterator;

/**
 *  The calculation of bond-pi Partial charge is calculated 
 *  determining the difference the Partial Pi Charge on atoms 
 *  A and B of a bond. Based in Gasteiger Charge.
 *  <p>This descriptor uses these parameters:
 * <table border="1">
 *   <tr>
 *     <td>Name</td>
 *     <td>Default</td>
 *     <td>Description</td>
 *   </tr>
 *   <tr>
 *     <td>bondPosition</td>
 *     <td>0</td>
 *     <td>The position of the target bond</td>
 *   </tr>
 * </table>
 *
 *
 * @author      Miguel Rojas
 * @cdk.created 2006-05-18
 * @cdk.module  qsarbond
 * @cdk.githash
 * @cdk.set     qsar-descriptors
 * @cdk.dictref qsar-descriptors:bondPartialPiCharge
 *
 * @see org.openscience.cdk.qsar.descriptors.atomic.PartialPiChargeDescriptor
 */
@TestClass(value="org.openscience.cdk.qsar.descriptors.bond.BondPartialPiChargeDescriptorTest")
public class BondPartialPiChargeDescriptor extends AbstractBondDescriptor {

    private GasteigerPEPEPartialCharges pepe = null;
    /**Number of maximum iterations*/
	private int maxIterations = -1;
    /**Number of maximum resonance structures*/
	private int maxResonStruc = -1;
	/** make a lone pair electron checker. Default true*/
	private boolean lpeChecker = true;

    private static final String[] descriptorNames = {"pepeB"};
    /**
     *  Constructor for the BondPartialPiChargeDescriptor object
     */
    public BondPartialPiChargeDescriptor() {  
    	pepe = new GasteigerPEPEPartialCharges();
    }


    /**
     *  Gets the specification attribute of the BondPartialPiChargeDescriptor
     *  object
     *
     *@return    The specification value
     */
    @TestMethod(value="testGetSpecification")
    public DescriptorSpecification getSpecification() {
        return new DescriptorSpecification(
            "http://www.blueobelisk.org/ontologies/chemoinformatics-algorithms/#bondPartialPiCharge",
            this.getClass().getName(),
            "$Id: 4a47581bde571789c1c0eae97275743af0f372a8 $",
            "The Chemistry Development Kit");
    }

    /**
     * This descriptor does have any parameter.
     */
    @TestMethod(value="testSetParameters_arrayObject")
    public void setParameters(Object[] params) throws CDKException {
        if (params.length > 3) 
            throw new CDKException("PartialPiChargeDescriptor only expects three parameter");
        
        if (!(params[0] instanceof Integer) )
                throw new CDKException("The parameter must be of type Integer");
	        maxIterations = (Integer) params[0];
	        
	    if(params.length > 1 && params[1] != null){
        	if (!(params[1] instanceof Boolean) )
                throw new CDKException("The parameter must be of type Boolean");
        	lpeChecker = (Boolean) params[1];
        }
	    
	    if(params.length > 2 && params[2] != null){
        	if (!(params[2] instanceof Integer) )
                throw new CDKException("The parameter must be of type Integer");
        	maxResonStruc = (Integer) params[2];
        }
    }


    /**
     *  Gets the parameters attribute of the BondPartialPiChargeDescriptor object.
     *
     *@return    The parameters value
     * @see #setParameters
     */
    @TestMethod(value="testGetParameters")
    public Object[] getParameters() {
        // return the parameters as used for the descriptor calculation
        Object[] params = new Object[3];
        params[0] = maxIterations;
        params[1] = lpeChecker;
        params[2] = maxResonStruc;
        return params;
    }

    @TestMethod(value="testNamesConsistency")
    public String[] getDescriptorNames() {
        return descriptorNames;
    }


    private DescriptorValue getDummyDescriptorValue(Exception e) {
        return new DescriptorValue(getSpecification(), getParameterNames(), getParameters(),
                new DoubleResult(Double.NaN), descriptorNames,e);
    }
    /**
     *  The method calculates the bond-pi Partial charge of a given bond
     *  It is needed to call the addExplicitHydrogensToSatisfyValency method from the class tools.HydrogenAdder.
     *
     *@param  ac                AtomContainer
     *@return                   return the sigma electronegativity
     */
    @TestMethod(value="testCalculate_IBond_IAtomContainer,testBondPiElectronegativityDescriptor,testBondPiElectronegativityDescriptor_Allyl_bromide")
    public DescriptorValue calculate(IBond bond, IAtomContainer ac)  {
    	// FIXME: for now I'll cache a few modified atomic properties, and restore them at the end of this method
    	Double originalCharge1 = bond.getAtom(0).getCharge();
    	String originalAtomtypeName1 = bond.getAtom(0).getAtomTypeName();
    	Integer originalNeighborCount1 = bond.getAtom(0).getFormalNeighbourCount();
    	IAtomType.Hybridization originalHybridization1 = bond.getAtom(0).getHybridization();
    	Integer originalValency1 = bond.getAtom(0).getValency();
    	Double originalCharge2 = bond.getAtom(1).getCharge();
    	String originalAtomtypeName2 = bond.getAtom(1).getAtomTypeName();
    	Integer originalNeighborCount2 = bond.getAtom(1).getFormalNeighbourCount();
    	IAtomType.Hybridization originalHybridization2 = bond.getAtom(1).getHybridization();
    	Integer originalValency2 = bond.getAtom(1).getValency();
        if (!isCachedAtomContainer(ac)) {
            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
                if (lpeChecker) {
                    LonePairElectronChecker lpcheck = new LonePairElectronChecker();
                    lpcheck.saturate(ac);
                }
            } catch (CDKException e) {
                return getDummyDescriptorValue(e);
            }

            if(maxIterations != -1)
    			pepe.setMaxGasteigerIters(maxIterations);
    		if(maxResonStruc != -1)
    			pepe.setMaxResoStruc(maxResonStruc);
	    	try {
	    		for (int i=0; i<ac.getAtomCount(); i++)
	    			ac.getAtom(i).setCharge(0.0);
	    		
	        	pepe.assignGasteigerPiPartialCharges(ac, true);
				for(Iterator<IBond> it = ac.bonds().iterator() ; it.hasNext(); ) {
					IBond bondi = it.next();
					double result = Math.abs(bondi.getAtom(0).getCharge()-bondi.getAtom(1).getCharge());
					cacheDescriptorValue(bondi, ac, new DoubleResult(result));
				}
	        } catch (Exception ex1) {
	            return getDummyDescriptorValue(ex1);
	        }
    	}
	    bond.getAtom(0).setCharge(originalCharge1);
	    bond.getAtom(0).setAtomTypeName(originalAtomtypeName1);
	    bond.getAtom(0).setHybridization(originalHybridization1);
	    bond.getAtom(0).setValency(originalValency1);
	    bond.getAtom(0).setFormalNeighbourCount(originalNeighborCount1);
	    bond.getAtom(1).setCharge(originalCharge2);
	    bond.getAtom(1).setAtomTypeName(originalAtomtypeName2);
	    bond.getAtom(1).setHybridization(originalHybridization2);
	    bond.getAtom(1).setValency(originalValency2);
	    bond.getAtom(1).setFormalNeighbourCount(originalNeighborCount2);
	    
        return getCachedDescriptorValue(bond) != null
                ? new DescriptorValue(getSpecification(), getParameterNames(), getParameters(),
                getCachedDescriptorValue(bond), descriptorNames)
                : null;
    }

	 /**
   * Gets the parameterNames attribute of the BondPartialPiChargeDescriptor object.
   *
   * @return    The parameterNames value
   */
    @TestMethod(value="testGetParameterNames")
    public String[] getParameterNames() {
    	String[] params = new String[3];
        params[0] = "maxIterations";
        params[1] = "lpeChecker";
        params[2] = "maxResonStruc";
        return params;
    }


  /**
   * Gets the parameterType attribute of the BondPartialPiChargeDescriptor object.
   *
   * @param  name  Description of the Parameter
   * @return       An Object of class equal to that of the parameter being requested
   */
    @TestMethod(value="testGetParameterType_String")
    public Object getParameterType(String name) {
    	if ("maxIterations".equals(name)) return Integer.MAX_VALUE;
    	if ("lpeChecker".equals(name)) return Boolean.TRUE;
    	if ("maxResonStruc".equals(name)) return Integer.MAX_VALUE;
        return null;
    }
}

