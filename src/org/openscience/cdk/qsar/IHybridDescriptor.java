package org.openscience.cdk.qsar;

import org.openscience.cdk.interfaces.IAtomContainer;

public interface IHybridDescriptor extends IDescriptor {
	
	    /** 
	     * Calculates the hybrid descriptor value for the given IAtom.
	     *
	     * @param  container        An {@link IAtomContainer} for which this descriptor should be
	     *                      calculated
	     * 
	     * @return              An object of {@link IAtomContainer} that contain the 
	     *                      calculated value as well as specification details
	     */
	    
	public enum Descriptor {
		REFRACTEDTOPOLOGICAL,
		REFRACTEDTOPOGRAPHIC,
		REFRACTEDTOPOGRAPHICNORM,
		ELECTROTOPOLOGICAL,
		ELECTROTOPOGRAPHIC,
		ELECTROTOPOGRAPHICNORM,
		LIPOTOPOLOGICAL,
		LIPOTOPOGRAPHIC,
		LIPOTOPOGRAPHICNORM,
		TOTALREFRACTEDTOPOLOGICAL,
		TOTALREFRACTEDTOPOGRAPHIC,
		TOTALELECTROTOPOLOGICAL,
		TOTALELECTROTOPOGRAPHIC,
		TOTALLIPOTOPOLOGICAL,
		TOTALLIPOTOPOGRAPHIC,
		REFRACTIINTRINSICVALUE,
		ELECTROINTRINSICVALUE,
		LIPOINTRINSICVALUE		
	}
	
	public DescriptorValue calculate(IAtomContainer atomcontainer);
	
	public IAtomContainer getAtomContainer();
	
	
	
	

}
