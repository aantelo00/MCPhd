package org.openscience.cdk.qsar.descriptors.hybrid.topographical;

import java.lang.reflect.Method;
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.qsar.DescriptorSpecification;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.DistanceToAtomDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.DoubleResultType;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.tools.AtomicProperties;
import org.openscience.cdk.tools.LoggingTool;
import org.openscience.cdk.tools.manipulator.BondManipulator;

public abstract class TopographicDescriptor implements IHybridDescriptor {
	
	protected LoggingTool logger;
	
	protected IAtomContainer atomContainer, container;
	protected IRingSet rs;
	protected String[] fragment; // estate fragments for each atom

	protected AtomicProperties ap; // needed to retrieve electronegativities

	public int[] alogpfrag; // alogp fragments for each atom (used to see which atoms have missing fragments)
	
	public TopographicDescriptor() throws CDKException {
		logger = new LoggingTool(this);

		try {
			ap = AtomicProperties.getInstance();
		} catch (Exception e) {
			logger.debug("Problem in accessing atomic properties. Can't calculate");
			throw new CDKException("Problem in accessing atomic properties. Can't calculate\n" + e.getMessage(), e);
		}
	}
	
	protected String UnassignedAtoms="";
	
	private void findUnassignedAtoms() {
		UnassignedAtoms="";

		for (int i = 0; i <= atomContainer.getAtomCount() - 1; i++) {
			if (alogpfrag[i]==0) UnassignedAtoms+=(i+1)+"("+fragment[i]+"),";
		}
	}

    protected void calculate(IAtomContainer atomContainer, String[] fragment, IRingSet rs) throws CDKException {
		this.atomContainer = atomContainer;
		this.fragment = fragment;
		this.rs = rs;

		alogpfrag = new int[atomContainer.getAtomCount()];

		for (int i = 0; i < atomContainer.getAtomCount(); i++) {

			alogpfrag[i] = 0;
			try {
				// instead of calling hardcoded methods here, use retrospection
				// and run all methods whos name start with 'calc' except for
				// 'calculate'. Nice :)
				@SuppressWarnings("unused")
				Method[] methods = this.getClass().getDeclaredMethods();
				if (fragment[i] instanceof String) {
					if(!atomContainer.getAtom(i).getSymbol().equals("H")){
					  	calcDescriptor(fragment[i], i);	
						@SuppressWarnings("unused")
					    double valorEstateIntAtomoOrigen = Double.valueOf(atomContainer.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE).toString());
					}
				}
				
			} catch (Exception e) {
				throw new CDKException(e.toString(), e);
			}
		} // end i atom loop

		logger.debug("\nFound fragments and frequencies ");
		this.findUnassignedAtoms();
	}
    
    public abstract void calcDescriptor(String cadena, int i);

    protected int getHAtomType(IAtom ai, List<IAtom> connectedAtoms)
	{
		//ai is the atom connected to a H atoms.
		//ai environment determines what is the H atom type
		//This procedure is applied only for carbons
		//i.e. H atom type 50 is never returned

		List<IAtom> ca;
		if (connectedAtoms == null)
			ca = atomContainer.getConnectedAtomsList(ai);
		else
			ca = connectedAtoms;

		// first check for alpha carbon:
		if (ai.getSymbol().equals("C") && !ai.getFlag(CDKConstants.ISAROMATIC)) {
			for (int j = 0; j <= ca.size() - 1; j++)
			{
				if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE && ((IAtom)ca.get(j)).getSymbol().equals("C")) { // single bonded
					List<IAtom> ca2 = atomContainer.getConnectedAtomsList((IAtom)ca.get(j));

					for (int k = 0; k <= ca2.size() - 1; k++)
					{
						IAtom ca2k = (IAtom)ca2.get(k);
						if (!ca2k.getSymbol().equals("C"))
						{
							if (atomContainer.getBond(((IAtom)ca.get(j)), ca2k).getOrder() != IBond.Order.SINGLE)
								return 51;

							if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)
									&& ca2k.getFlag(CDKConstants.ISAROMATIC)) {
								if (inSameAromaticRing(atomContainer, ((IAtom)ca.get(j)), ca2k,	rs))
								{
									return 51;
								}
							}
						} // end !ca2[k].getSymbol().equals("C"))
					} // end k loop
				} // end if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
			}// end j loop
		} // end if(ai.getSymbol().equals("C") && !ai.getFlag(CDKConstants.ISAROMATIC))

		List<IBond> bonds = atomContainer.getConnectedBondsList(ai);
		int doublebondcount = 0;
		int triplebondcount = 0;
		String hybrid = "";

		for (int j = 0; j < bonds.size(); j++)
		{
			if (((IBond)bonds.get(j)).getOrder() == IBond.Order.DOUBLE)
				doublebondcount++;
			else
				if (((IBond)bonds.get(j)).getOrder() == IBond.Order.TRIPLE)
					triplebondcount++;
		}

		if (doublebondcount == 0 && triplebondcount == 0)
			hybrid = "sp3";
		else
			if (doublebondcount == 1 && triplebondcount == 0)
				hybrid = "sp2";
			else
				if (doublebondcount == 2 || triplebondcount == 1)
					hybrid = "sp";
		int OxNum = 0;
		int XCount = 0;

		for (int j = 0; j < ca.size(); j++)
		{
			//String s = ((IAtom)ca.get(j)).getSymbol();
			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))
			if (ap.getNormalizedElectronegativity(((IAtom)ca.get(j)).getSymbol()) > 1)
			{
				List<IBond> bonds2 = atomContainer.getConnectedBondsList(((IAtom)ca.get(j)));
				boolean HaveDouble = false;
				for (int k = 0; k < bonds2.size(); k++)
				{
					if (((IBond)bonds2.get(k)).getOrder() == IBond.Order.DOUBLE)
					{
						HaveDouble = true;
						break;
					}
				}
				if (HaveDouble && ((IAtom)ca.get(j)).getSymbol().equals("N"))
					OxNum += 2; // C-N bond order for pyridine type N's is considered to be 2
					else
						OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder());
			}
			List<IAtom> ca2 = atomContainer.getConnectedAtomsList(((IAtom)ca.get(j)));

			for (int k = 0; k < ca2.size(); k++)
			{
				String s2 = ((IAtom)ca2.get(k)).getSymbol();
				if (!s2.equals("C"))
					XCount++;
			}
		}// end j loop

		if (OxNum == 0)
		{
			if (hybrid.equals("sp3"))
			{
				if (XCount == 0)
					return 46;
				else if (XCount == 1)
					return 52;
				else if (XCount == 2)
					return 53;
				else if (XCount == 3)
					return 54;
				else if (XCount >= 4)
					return 55;
			}
			else if (hybrid.equals("sp2"))
				return 47;
		}
		else if (OxNum == 1 && hybrid.equals("sp3"))
			return 47;
		else if ((OxNum == 2 && hybrid.equals("sp3"))
				|| (OxNum == 1 && hybrid.equals("sp2"))
				|| (OxNum == 0 && hybrid.equals("sp")))
			return 48;
		else if ((OxNum == 3 && hybrid.equals("sp3"))
				|| (OxNum >= 2 && hybrid.equals("sp2"))
				|| (OxNum >= 1 && hybrid.equals("sp")))
			return 49;

		return(0);
	}
	
	protected boolean inSameAromaticRing(IAtomContainer atomContainer, IAtom atom1,
			IAtom atom2, IRingSet rs) {
		boolean SameRing = false;

		for (int i = 0; i <= rs.getAtomContainerCount() - 1; i++) {
			IRing r = (IRing)rs.getAtomContainer(i);

			if (!r.getFlag(CDKConstants.ISAROMATIC))
				continue;

			// ArrayList al=new ArrayList();

			boolean HaveOne = false;
			boolean HaveTwo = false;

			for (int j = 0; j <= r.getAtomCount() - 1; j++) {
				if (atomContainer.getAtomNumber(r.getAtom(j)) == atomContainer.getAtomNumber(atom1))
					HaveOne = true;
				if (atomContainer.getAtomNumber(r.getAtom(j)) == atomContainer.getAtomNumber(atom2))
					HaveTwo = true;
			}

			if (HaveOne && HaveTwo) {
				SameRing = true;
				return SameRing;
			}

		} // end ring for loop

		return SameRing;
	}

	/**
	 * The AlogP descriptor.
	 *
	 * TODO Ideally we should explicit H addition should be cached
	 *
	 * @param atomContainer the molecule to calculate on
	 * @return the result of the calculation
	 */
	
	@TestMethod("testCalculate_IAtomContainer,testChloroButane")
	public abstract DescriptorValue calculate(IAtomContainer atomContainer);

	protected void selectAromaticRing(IRingSet rs2) {

		boolean aromatic;
		for(int i = 0; i < rs2.getAtomContainerCount(); i++){
			aromatic = true;
			IAtomContainer contAux = rs2.getAtomContainer(i);
			for(IBond bond : contAux.bonds()){
				if(!bond.getFlag(CDKConstants.ISAROMATIC)){
					aromatic = false;
					break;
				}    				
			}
			rs2.getAtomContainer(i).setFlag(CDKConstants.ISAROMATIC, aromatic);
		}

	}

	protected double updateRefractedTopographicValues() throws CDKException {

		IAtomicDescriptor distance = new DistanceToAtomDescriptor();
		double efecto = 0.0, total = 0.0;
		//Object[] arrAux = new Object[1];

		int nAtomos = container.getAtomCount();
		for (int i = 0; i < nAtomos; i++) {
			efecto = 0.0;
			IAtom atomoOrigen = container.getAtom(i);
			if(!atomoOrigen.getSymbol().equals("H")){
				double valorEstateIntAtomoOrigen = Double.valueOf(atomoOrigen.getProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE).toString());
				for(int j = 0; j < nAtomos; j++){
					IAtom atomoDestino = container.getAtom(j);
					if(i != j && !atomoDestino.getSymbol().equals("H")){

						//arrAux[0] = j;
						Object[] objs = {j};
						distance.setParameters(objs);
						//(IntegerArrayResult)descriptor.calculate(containersList.get(0)).getValue();
						DescriptorValue value =  distance.calculate(atomoOrigen, container);
						double dist = ((DoubleResult) value.getValue()).doubleValue();
						if(dist != 0){
							double valorEstateIntAtomoDestino = Double.valueOf(atomoDestino.getProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE).toString());
							efecto += (valorEstateIntAtomoOrigen - valorEstateIntAtomoDestino)/Math.pow(dist, 2);

						}
					}

				}

				atomoOrigen.setProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC, valorEstateIntAtomoOrigen + efecto);
				total += valorEstateIntAtomoOrigen + efecto;

			}
		}
		return total;
	}

	protected DescriptorValue getDummyDescriptorValue(Exception e) {
		DoubleArrayResult results = new DoubleArrayResult();
		results.add(Double.NaN);
		results.add(Double.NaN);
		results.add(Double.NaN);
		return new DescriptorValue(getSpecification(), getParameterNames(),
				getParameters(), results, getDescriptorNames(), e);
	}

	/**
	 * Returns the specific type of the DescriptorResult object.
	 * <p/>
	 * The return value from this method really indicates what type of result will
	 * be obtained from the {@link org.openscience.cdk.qsar.DescriptorValue} object. Note that the same result
	 * can be achieved by interrogating the {@link org.openscience.cdk.qsar.DescriptorValue} object; this method
	 * allows you to do the same thing, without actually calculating the descriptor.
	 *
	 * @return an object that implements the {@link org.openscience.cdk.qsar.result.IDescriptorResult} interface indicating
	 *         the actual type of values returned by the descriptor in the {@link org.openscience.cdk.qsar.DescriptorValue} object
	 */
	
	@TestMethod("testGetDescriptorResultType")
	protected IDescriptorResult getDescriptorResultType() {
		return new DoubleResultType();
	}

	//Ver como adecuar esta especificaci�n.
	@TestMethod("testGetSpecification")
	public DescriptorSpecification getSpecification() {
		return new DescriptorSpecification(
				"http://www.blueobelisk.org/ontologies/chemoinformatics-algorithms/#ALOGP",
				this.getClass().getName(),
				"$Id$",
		"The Chemistry Development Kit");
	}

	@TestMethod("testGetParameterNames")
	public String[] getParameterNames() {
		return new String[0];
	}

	@TestMethod("testGetParameterType_String")
	public Object getParameterType(String name) {
		return null;
	}

	@TestMethod("testSetParameters_arrayObject")
	public void setParameters(Object[] params) throws CDKException {
	}

	@TestMethod("testGetParameters")
	public Object[] getParameters() {
		return null;
	}

	@TestMethod(value="testNamesConsistency")
	public String[] getDescriptorNames() {
		return new String[]{"Refracto topographic descriptor"};
	}

	public IAtomContainer getAtomContainer(){
		return container;
	}
	
	//public abstract void calcAtomC(int i);
	//protected abstract void calcGroup001_005(int i);
	//protected abstract void calcGroup002_006_007(int i);
	//protected abstract void calcGroup003_008_009_010(int i);
	//protected abstract void calcGroup004_011_to_014(int i);
	//protected abstract void calcGroup015(int i);
	//protected abstract void calcGroup016_018_036_037(int i);
	//protected abstract void calcGroup017_019_020_038_to_041(int i);
	//protected abstract void calcGroup021_to_023_040(int i);
	//protected abstract void calcGroup024_027_030_033_042(int i);
	//protected abstract void calcGroup025_026_028_029_031_032_034_035_043_044(int i);
	//protected abstract void calcGroup056_57(int i);
	//protected abstract void calcGroup058_61(int i);
	//protected abstract void calcGroup059_060_063(int i);
	//protected abstract void calcGroup066_to_079(int i);
	//protected abstract void calcGroup081_to_085(int i);
	//protected abstract void calcGroup086_to_090(int i);
	//protected abstract void calcGroup091_to_095(int i);
	//protected abstract void calcGroup096_to_100(int i);
	//protected abstract void calcGroup101_to_104(int i);
	//protected abstract void calcGroup106(int i);
	//protected abstract void calcGroup107(int i);
	//protected abstract void calcGroup108(int i);
	//protected abstract void calcGroup109(int i);
	//protected abstract void calcGroup110(int i);
	//protected abstract void calcGroup111(int i);
	//protected abstract void calcGroup116_117_120(int i);
    //protected abstract void calcGroup118_119(int i);

}
