package org.openscience.cdk.qsar.descriptors.hybrid.topographical;

import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.atomtype.EStateAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.DistanceToAtomDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.BondManipulator;

public class RefractedTopographicDescriptor extends TopographicDescriptor {

	final static double[] refracval = new double[121]; // coefficients for refractivity model
	static {
		// fragments for AMR from Viswanadhan et al., 1989
		refracval[1]=2.968;
		refracval[2]=2.9116;
		refracval[3]=2.8028;
		refracval[4]=2.6205;
		refracval[5]=3.015;
		refracval[6]=2.9244;
		refracval[7]=2.6329;
		refracval[8]=2.504;
		refracval[9]=2.377;
		refracval[10]=2.5559;
		refracval[11]=2.303;
		refracval[12]=2.3006;
		refracval[13]=2.9627;
		refracval[14]=2.3038;
		refracval[15]=3.2001;
		refracval[16]=4.2654;
		refracval[17]=3.9392;
		refracval[18]=3.6005;
		refracval[19]=4.487;
		refracval[20]=3.2001;
		refracval[21]=3.4825;
		refracval[22]=4.2817;
		refracval[23]=3.9556;
		refracval[24]=3.4491;
		refracval[25]=3.8821;
		refracval[26]=3.7593;
		refracval[27]=2.5009;
		refracval[28]=2.5;
		refracval[29]=3.0627;
		refracval[30]=2.5009;
		refracval[31]=0;
		refracval[32]=2.6632;
		refracval[33]=3.4671;
		refracval[34]=3.6842;
		refracval[35]=2.9372;
		refracval[36]=4.019;
		refracval[37]=4.777;
		refracval[38]=3.9031;
		refracval[39]=3.9964;
		refracval[40]=3.4986;
		refracval[41]=3.4997;
		refracval[42]=2.7784;
		refracval[43]=2.6267;
		refracval[44]=2.5;
		refracval[45]=0;
		refracval[46]=0.8447;
		refracval[47]=0.8939;
		refracval[48]=0.8005;
		refracval[49]=0.832;
		refracval[50]=0.8;
		refracval[51]=0.8188;
		refracval[52]=0.9215;
		refracval[53]=0.9769;
		refracval[54]=0.7701;
		refracval[55]=0;
		refracval[56]=1.7646;
		refracval[57]=1.4778;
		refracval[58]=1.4429;
		refracval[59]=1.6191;
		refracval[60]=1.3502;
		refracval[61]=1.945;
		refracval[62]=0;
		refracval[63]=0;
		refracval[64]=11.1366;
		refracval[65]=13.1149;
		refracval[66]=2.6221;
		refracval[67]=2.5;
		refracval[68]=2.898;
		refracval[69]=3.6841;
		refracval[70]=4.2808;
		refracval[71]=3.6189;
		refracval[72]=2.5;
		refracval[73]=2.7956;
		refracval[74]=2.7;
		refracval[75]=4.2063;
		refracval[76]=4.0184;
		refracval[77]=3.0009;
		refracval[78]=4.7142;
		refracval[79]=0;
		refracval[80]=0;
		refracval[81]=0.8725;
		refracval[82]=1.1837;
		refracval[83]=1.1573;
		refracval[84]=0.8001;
		refracval[85]=1.5013;
		refracval[86]=5.6156;
		refracval[87]=6.1022;
		refracval[88]=5.9921;
		refracval[89]=5.3885;
		refracval[90]=6.1363;
		refracval[91]=8.5991;
		refracval[92]=8.9188;
		refracval[93]=8.8006;
		refracval[94]=8.2065;
		refracval[95]=8.7352;
		refracval[96]=13.9462;
		refracval[97]=14.0792;
		refracval[98]=14.073;
		refracval[99]=12.9918;
		refracval[100]=13.3408;
		refracval[101]=0;
		refracval[102]=0;
		refracval[103]=0;
		refracval[104]=0;
		refracval[105]=0;
		refracval[106]=7.8916;
		refracval[107]=7.7935;
		refracval[108]=9.4338;
		refracval[109]=7.7223;
		refracval[110]=5.7558;
		refracval[111]=0;
		refracval[112]=0;
		refracval[113]=0;
		refracval[114]=0;
		refracval[115]=0;
		refracval[116]=5.5306;
		refracval[117]=5.5152;
		refracval[118]=6.836;
		refracval[119]=10.0101;
		refracval[120]=5.2806;
	}

	
	public RefractedTopographicDescriptor() throws CDKException {
		super();
	}
	
	@Override
	public void calcDescriptor(String cadena, int i){
		if(cadena.equals("sCH3") || cadena.equals("ssCH2") || cadena.equals("sssCH") || cadena.equals("ssssC") ||
		   cadena.equals("dCH2") || cadena.equals("dsCH") || cadena.equals("dssC") || cadena.equals("tCH")	||
		   cadena.equals("tsC") || cadena.equals("aaCH") || cadena.equals("saaC") || cadena.equals("aaaC") ||
		   cadena.equals("daaC") || cadena.equals("ddC"))
			calcAtomC(i);
		else if(cadena.equals("sOH"))
			calcGroup056_57(i);
		else if(cadena.equals("sOm") || cadena.equals("dO"))
			calcGroup058_61(i);
		else if(cadena.equals("ssO") || cadena.equals("aaO") || cadena.equals("aaOp"))
			calcGroup059_060_063(i);
		else if(cadena.equals("sNH2") || cadena.equals("aaNH") || cadena.equals("saaN") || cadena.equals("ssNH") ||
				cadena.equals("sssN") || cadena.equals("aaN") || cadena.equals("aaaN") || cadena.equals("aaNm") ||
				cadena.equals("ssdNp") || cadena.equals("dssNp") || cadena.equals("tN") || cadena.equals("dNH") || 
				cadena.equals("dsN"))
			calcGroup066_to_079(i);
		else if(cadena.equals("sF"))
			calcGroup081_to_085(i);
		else if(cadena.equals("sCl"))
			calcGroup086_to_090(i);
		else if(cadena.equals("sBr"))
			calcGroup091_to_095(i);
		else if(cadena.equals("sI"))
			calcGroup096_to_100(i);
		else if(cadena.equals("F") || cadena.equals("Cl") || cadena.equals("Br") || cadena.equals("I"))
			calcGroup101_to_104(i);
		else if(cadena.equals("sSH"))
			calcGroup106(i);
		else if(cadena.equals("ssS") || cadena.equals("aaS") || cadena.equals("sSm"))
			calcGroup107(i);
		else if(cadena.equals("dS"))
			calcGroup108(i);
		else if(cadena.equals("dssS"))
			calcGroup109(i);
		else if(cadena.equals("ddssS"))
			calcGroup110(i);
		else if(cadena.equals("ssssSi"))
			calcGroup111(i);
		else if(cadena.equals("dsssP"))
			calcGroup116_117_120(i);
		else if(cadena.equals("sssP"))
			calcGroup118_119(i);
	}
	
	public void calcAtomC(int i) {
		int htype = 0;
		int atomos5 = 0;
		int CarbonCount = 0;
		int HeteroCount = 0;
		double perturbacion = 0.0;
		boolean HaveCdX = false;
		boolean HaveCsX = false;
		boolean HaveCsAr = false;
		boolean HaveCtX = false;
		IAtom atomoPrincipal = atomContainer.getAtom(i);
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomoPrincipal);
		
	    String cadena =  fragment[i];
	    
		switch (cadena) {
		case "sCH3": {
			perturbacion = 0.0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j < ca.size(); j++) {
				IAtom atomoAux = ca.get(j);
				if (atomoAux.getSymbol().equals("C"))
					alogpfrag[i] = 1;
				/*
				 * IAtom atomoAux = ca.get(j); if (atomoAux.equals("C"))
				 * alogpfrag[i] = 1;
				 */
				else if (atomoAux.getSymbol().equals("H"))
					perturbacion += refracval[htype];
				else {
					alogpfrag[i] = 5;
					atomos5++;
				}
			}
			if (atomos5 == 0)
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[1] + perturbacion);
			else
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[5] + perturbacion);
			break;
		}
		case "ssCH2": {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j < ca.size(); j++)
				if (((IAtom) ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else if (((IAtom) ca.get(j)).getSymbol().equals("H"))
					perturbacion += refracval[htype];
				else
					HeteroCount++;
			if (CarbonCount == 2 && HeteroCount == 0) {
				alogpfrag[i] = 2;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[2] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 1) {
				alogpfrag[i] = 6;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[6] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 2) {
				alogpfrag[i] = 7;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[7] + perturbacion);
			}
			break;
		}
		case "sssCH": {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j <= ca.size() - 1; j++)
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					perturbacion += refracval[htype];
			    else
					HeteroCount++;
			if (CarbonCount == 3 && HeteroCount == 0) {
				alogpfrag[i] = 3;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[3] + perturbacion);
			} else if (CarbonCount == 2 && HeteroCount == 1) {
				alogpfrag[i] = 8;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[8] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 2) {
				alogpfrag[i] = 9;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[9] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 3) {
				alogpfrag[i] = 10;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[10] + perturbacion);
			}
			break;
		}
		case "ssssC": {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			for (int j = 0; j <= ca.size() - 1; j++)
				if (((IAtom) ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else
					HeteroCount++;
			if (CarbonCount == 4 && HeteroCount == 0) {
				alogpfrag[i] = 4;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[4] + perturbacion);
			} else if (CarbonCount == 3 && HeteroCount == 1) {
				alogpfrag[i] = 11;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[11] + perturbacion);
			} else if (CarbonCount == 2 && HeteroCount == 2) {
				alogpfrag[i] = 12;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[12] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 3) {
				alogpfrag[i] = 13;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[13] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 4) {
				alogpfrag[i] = 14;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[14] + perturbacion);
			}
			break;
		}
		case "dCH2": {
			perturbacion = 0.0;
			htype = getHAtomType(atomContainer.getAtom(i), null);
			alogpfrag[i] = 15;
			perturbacion += 2 * refracval[htype];
			atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[15] + perturbacion);
			break;
		}
		case "dsCH": {
			perturbacion = 0.0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			perturbacion += refracval[htype];
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom) ca.get(j)).getSymbol().equals("H")) {
					continue;
				}
				if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j)))
						.getOrder() == IBond.Order.SINGLE) {
					if (!((IAtom) ca.get(j)).getSymbol().equals("C"))
						HaveCsX = true;
					if (((IAtom) ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
						HaveCsAr = true;
				} else if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j))).getOrder() == IBond.Order.DOUBLE) 
					if (!((IAtom) ca.get(j)).getSymbol().equals("C"))
						HaveCdX = true;
			}
			if (HaveCdX) {
				if (HaveCsAr) {
					alogpfrag[i] = 37;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[37] + perturbacion);
				} else {
					alogpfrag[i] = 36;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[36] + perturbacion);
				}
			} else {
				if (HaveCsX) {
					alogpfrag[i] = 18;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[18] + perturbacion);
				} else {
					alogpfrag[i] = 16;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[16] + perturbacion);
				}
			}
			break;
		}
		case "dssC": {
			HaveCdX = false;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j))).getOrder() == IBond.Order.SINGLE) {
					if (((IAtom) ca.get(j)).getSymbol().equals("C")) {
					} else {
					}
					if (!((IAtom) ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
					} else {
					}
				} else if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j))).getOrder() == IBond.Order.DOUBLE)
					if (!((IAtom) ca.get(j)).getSymbol().equals("C"))
						HaveCdX = true;
			}
			if (HaveCdX) {
				if (HaveCsAr) {
					alogpfrag[i] = 37;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[37] + perturbacion);
				} else {
					alogpfrag[i] = 36;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[36] + perturbacion);
				}
			} else {
				if (HaveCsX) {
					alogpfrag[i] = 18;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[18] + perturbacion);
				} else {
					alogpfrag[i] = 16;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[16] + perturbacion);
				}
			}
			break;
		}
		case "tCH": {
			alogpfrag[i] = 21;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[21] + refracval[htype]);
			break;
		}
		case "ddC": {
			if (((IAtom) ca.get(0)).getSymbol().equals("C")
					&& ((IAtom) ca.get(1)).getSymbol().equals("C")) {
				alogpfrag[i] = 22;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE,refracval[22]);
			} else if (!((IAtom) ca.get(0)).getSymbol().equals("C")
					&& !((IAtom) ca.get(1)).getSymbol().equals("C")) {
				alogpfrag[i] = 40;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[40]);
			}
			break;
		}
		case "tsC": {
			HaveCtX = false;
			HaveCsX = false;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j))).getOrder() == IBond.Order.SINGLE)
					if (!((IAtom) ca.get(j)).getSymbol().equals("C"))
						HaveCsX = true;
					else if (atomContainer.getBond(atomoPrincipal, ((IAtom) ca.get(j))).getOrder() == IBond.Order.TRIPLE)
						if (!((IAtom) ca.get(j)).getSymbol().equals("C"))
							HaveCtX = true;
			}
			if (HaveCtX && !HaveCsX) {
				alogpfrag[i] = 40;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[40]);
			} else if (HaveCsX) {// #C-X
				alogpfrag[i] = 23;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[23]);
			} else if (!HaveCsX) { // #C-R
				alogpfrag[i] = 22;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[22]);
			}
			break;
		}
		case "aaCH": {
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			IAtom ca0;
			IAtom ca1;
			if (((IAtom) ca.get(0)).getSymbol().equals("H")) {
				ca0 = (IAtom) ca.get(1);
				ca1 = (IAtom) ca.get(2);
			} else {
				if (((IAtom) ca.get(1)).getSymbol().equals("H")) {
					ca0 = (IAtom) ca.get(0);
					ca1 = (IAtom) ca.get(2);
				} else {
					ca0 = (IAtom) ca.get(0);
					ca1 = (IAtom) ca.get(1);
				}
			}
			if (ca0.getSymbol().equals("C") && ca1.getSymbol().equals("C")) {
				alogpfrag[i] = 24;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[24] + refracval[htype]);
				break;
			}
			List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
			boolean HaveDouble1 = false;
			for (int k = 0; k <= bonds.size() - 1; k++)
				if (((IBond) bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble1 = true;
					break;
				}
			bonds = atomContainer.getConnectedBondsList(ca1);
			boolean HaveDouble2 = false;
			for (int k = 0; k <= bonds.size() - 1; k++)
				if (((IBond) bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble2 = true;
					break;
				}
			if (!(ca0).getSymbol().equals("C") && !((IAtom) ca.get(1)).getSymbol().equals("C")) {
				if (HaveDouble1 && HaveDouble2) { // X--CH--X
					alogpfrag[i] = 30;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[30] + refracval[htype]);
				} else { // X--CH...X
					alogpfrag[i] = 42;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[42] + refracval[htype]);
				}

			} else if (ca0.getSymbol().equals("C")&& !ca1.getSymbol().equals("C") || (!ca0.getSymbol().equals("C") && ca1.getSymbol().equals("C"))) {
				if (HaveDouble1 && HaveDouble2) { // R--CH--X
					alogpfrag[i] = 27;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[27] + refracval[htype]);
				} else {// R--CH...X
					alogpfrag[i] = 33;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE,refracval[33] + refracval[htype]);
				}
			}
		}
		}
	    // quitar daaC ***********Importante**************
		if (fragment[i].equals("saaC") || fragment[i].equals("aaaC") || fragment[i].equals("daaC")) {
			IAtom[] sameringatoms = new IAtom[2];
			IAtom nonringatom = atomContainer.getBuilder().newInstance(IAtom.class);
			int sameringatomscount = 0;
			for (int j = 0; j <= ca.size() - 1; j++)
				if (inSameAromaticRing(atomContainer, atomoPrincipal, ((IAtom)ca.get(j)), rs))
					sameringatomscount++;
			if (sameringatomscount == 2) {
				int count = 0;
				for (int j = 0; j <= ca.size() - 1; j++) {
					if (inSameAromaticRing(atomContainer, atomoPrincipal, ((IAtom)ca.get(j)), rs)) {
						sameringatoms[count] = (IAtom)ca.get(j);
						count++;
					} else 
						nonringatom = (IAtom)ca.get(j);
				}
			} else {
				sameringatoms[0] = (IAtom)ca.get(0);
			    sameringatoms[1] = (IAtom)ca.get(1);
			    nonringatom = (IAtom)ca.get(2);
			}

			List<IBond> bonds = atomContainer.getConnectedBondsList(sameringatoms[0]);
			boolean HaveDouble1 = false;
			for (int k = 0; k <= bonds.size() - 1; k++)
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble1 = true;
					break;
				}
			bonds = atomContainer.getConnectedBondsList(sameringatoms[1]);
			boolean HaveDouble2 = false;
			for (int k = 0; k <= bonds.size() - 1; k++)
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble2 = true;
					break;
				}
			if (!sameringatoms[0].getSymbol().equals("C") && !sameringatoms[1].getSymbol().equals("C")) {
				if (HaveDouble1 && HaveDouble2) { // X--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						alogpfrag[i] = 31;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[31]);
					} else { // X--CX--X
						alogpfrag[i] = 32;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[32]);
					}
				} else {
					if (nonringatom.getSymbol().equals("C")) { // X--CR..X
						alogpfrag[i] = 43;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[43]);
					} else { // X--CX...X
						alogpfrag[i] = 44;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[44]);
					}
				}
			} else if (sameringatoms[0].getSymbol().equals("C") && sameringatoms[1].getSymbol().equals("C")) {
				if (nonringatom.getSymbol().equals("C")) {// R--CR--R
					alogpfrag[i] = 25;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[25]);
				} else { // R--CX--R
					alogpfrag[i] = 26;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[26]);
				}
			} else if ((sameringatoms[0].getSymbol().equals("C") && !sameringatoms[1].getSymbol().equals("C"))
			           || (!sameringatoms[0].getSymbol().equals("C") && sameringatoms[1].getSymbol().equals("C"))) {
				if (HaveDouble1 && HaveDouble2) { // R--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						alogpfrag[i] = 28;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[28]);
					} else { // R--CX--X
						alogpfrag[i] = 29;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[29]);
					}
				} else {
					if (nonringatom.getSymbol().equals("C")) { // R--CR..X
						alogpfrag[i] = 34;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[34]);
					} else { // R--CX...X
						alogpfrag[i] = 35;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[35]);
					}
				}
		    }
		}
	}
	
	public void calcGroup056_57(int i) {
		// 56: O in =O
		// 57: O in phenol, enol, and carboxyl
		// enol : compound containing a hydroxyl group bonded to a carbon atom
		// that in turn forms a double bond with another carbon atom.
		// enol = HO-C=C-
		// carboxyl= HO-C(=O)-

		if (!fragment[i].equals("sOH"))
			return;
		
		java.util.List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		if (ca0.getSymbol().equals("H"))
			ca0 = (IAtom)ca.get(1);
		if (ca0.getFlag(CDKConstants.ISAROMATIC)) { // phenol
			alogpfrag[i] = 57;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[57] + refracval[50]);
			return;
		}
		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		for (int j = 0; j <= ca2.size() - 1; j++) {
			if (atomContainer.getBond((IAtom)ca2.get(j), ca0).getOrder() == IBond.Order.DOUBLE) {
				alogpfrag[i] = 57;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[57] + refracval[50]);
				return;
			}
		}
		alogpfrag[i] = 56;
		atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[56] + refracval[50]);
	}

	public void calcGroup058_61(int i) {
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		if (fragment[i].equals("sOm")) {
			if (ca0.getSymbol().equals("N") && ca0.getFormalCharge() == 1) {
				alogpfrag[i] = 61;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[61]);
			} else {
				alogpfrag[i] = 62;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[62]);
			}
		} else if (fragment[i].equals("dO")) {
			if (ca0.getSymbol().equals("N") && ca0.getFormalCharge() == 1) {
				alogpfrag[i] = 61;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[61]);
			} else {
				alogpfrag[i] = 58;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[58]);
			}
		}
	}

	public void calcGroup059_060_063(int i) {
		if (!fragment[i].equals("ssO") && !fragment[i].equals("aaO") && !fragment[i].equals("aaOp"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		IAtom ca1 = (IAtom)ca.get(1);
		if (fragment[i].equals("ssO")) {
			if (ca0.getFlag(CDKConstants.ISAROMATIC) || ca1.getFlag(CDKConstants.ISAROMATIC)) {
				alogpfrag[i] = 60;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[60]);
			} else {
				for (int j = 0; j <= ca.size() - 1; j++) {
					List<IAtom> ca2 = atomContainer.getConnectedAtomsList(((IAtom)ca.get(j)));
					for (int k = 0; k <= ca2.size() - 1; k++) {
						if (atomContainer.getBond(((IAtom)ca.get(j)), (IAtom)ca2.get(k)).getOrder() == IBond.Order.DOUBLE) {
							if (!((IAtom)ca2.get(k)).getSymbol().equals("C")) {
								atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[60]);
								return;
							}
						}
					}
				}
				if (ca0.getSymbol().equals("O") || ca1.getSymbol().equals("O")) {
					alogpfrag[i] = 63;
					atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[63]);
				} else {
					alogpfrag[i] = 59;
					atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[59]);
				}
			}
		
			//Quitar aaOp ***********importante*********	
		} else if (fragment[i].equals("aaO") || fragment[i].equals("aaOp")) {
			alogpfrag[i] = 60;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[60]);
		}
	}

	public void calcGroup066_to_079(int i)
	{
		int NAr = 0;
		int NAl = 0;
		double perturbacion = 0.0;
		IAtom ai = atomContainer.getAtom(i);
		if (!ai.getSymbol().equals("N"))
			return;
		
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		for (int j = 0; j <= ca.size() - 1; j++) {
			if (((IAtom)ca.get(j)).getSymbol().equals("H"))
				continue;
			if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
				NAr++;
			else
				NAl++;
		}
		for (int j = 0; j <= ca.size() - 1; j++){
			if (((IAtom)ca.get(j)).getSymbol().equals("H"))
				continue;
			List<IAtom> ca2 = atomContainer.getConnectedAtomsList((IAtom)ca.get(j));
			for (int k = 0; k <= ca2.size() - 1; k++) {
				IAtom ca2k = (IAtom)ca2.get(k);
				if (atomContainer.getAtomNumber(ca2k) != i) {
					if (!ca2k.getSymbol().equals("C")) {
						if (!ca2k.getFlag(CDKConstants.ISAROMATIC)
								&& !((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)
								&& !ai.getFlag(CDKConstants.ISAROMATIC)) {
							if (atomContainer.getBond(((IAtom)ca.get(j)), ca2k).getOrder() == IBond.Order.DOUBLE) {
								alogpfrag[i] = 72;
								ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[72]);
								return;
							}
						}
					}
				}
			}
		}
		if (fragment[i].equals("sNH2")){
			IAtom ca0 = null;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					continue;
				else {
					ca0 = (IAtom)ca.get(j);
					break;
				}
			}
			if (ca0.getFlag(CDKConstants.ISAROMATIC) || !ca0.getSymbol().equals("C")) {
				alogpfrag[i] = 69;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[69] + 2 * refracval[50]);
			}
			else {
				alogpfrag[i] = 66;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[66] + 2 * refracval[50]);
			}
		}
		else if (fragment[i].equals("aaNH") || fragment[i].equals("saaN")) { // R...NH...R
			perturbacion = 0.0;
			alogpfrag[i] = 73;
			if (fragment[i].equals("aaNH")) {
				perturbacion = refracval[50];
			}
			ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[73] + perturbacion);
		}
		else if (fragment[i].equals("ssNH")) {
			if (NAr == 2 && NAl == 0) { // Ar2NH
				alogpfrag[i] = 73;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[73] + refracval[50]);
			} else if (NAr == 1 && NAl == 1) { // Ar-NH-Al
				alogpfrag[i] = 70;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[70] + refracval[50]);
			} else if (NAr == 0 && NAl == 2) { // Al2NH
				alogpfrag[i] = 67;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[67] + refracval[50]);
			}
		}
		else if (fragment[i].equals("sssN")) {
			if ((NAr == 3 && NAl == 0) || (NAr == 2 && NAl == 1)) { // Ar3N & // Ar2NAl
				alogpfrag[i] = 73;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[73]);
			} else if (NAr == 1 && NAl == 2) {
				alogpfrag[i] = 71;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[71]);
			} else if (NAr == 0 && NAl == 3) {
				alogpfrag[i] = 68;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[68]);
			}
		}
		
		//arreglar esto, quitar aaaN y aaNm  *****importante********
		else if (fragment[i].equals("aaN") || fragment[i].equals("aaaN") || fragment[i].equals("aaNm")) {
			alogpfrag[i] = 75;
			ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[75]);
		} else if (fragment[i].equals("ssdNp") || fragment[i].equals("dssNp")) {
			boolean HaveSsOm = false;
			boolean HaveSdO = false;
			boolean Ar = false;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (fragment[atomContainer.getAtomNumber(((IAtom)ca.get(j)))].equals("sOm")) {
					HaveSsOm = true;
				} else if (fragment[atomContainer.getAtomNumber(((IAtom)ca.get(j)))].equals("dO")) {
					HaveSdO = true;
				} else {
					if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
						Ar = true;
					}
				}
			}
			if (HaveSsOm && HaveSdO && Ar) {
				alogpfrag[i] = 76;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[76]);
			} else if (HaveSsOm && HaveSdO && !Ar) {
				alogpfrag[i] = 77;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[77]);
			} else {
				alogpfrag[i] = 79;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[79]);
			}

		} else if (fragment[i].equals("tN")) {
			IAtom ca0 = (IAtom)ca.get(0);
			//Quitar  !ca0.getSymbol().equals("C") ***********Importante**************
			if (ca0.getSymbol().equals("C") || !ca0.getSymbol().equals("C")) { // R#N
				alogpfrag[i] = 74;
				ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[74]);
			}
		} else if (fragment[i].equals("dNH") || fragment[i].equals("dsN")) { // test for RO-NO
			if (fragment[i].equals("dsN")) {
				IAtom ca0 = (IAtom)ca.get(0);
				IAtom ca1 = (IAtom)ca.get(1);
				if (ca0.getSymbol().equals("O")
						&& ca1.getSymbol().equals("O")) {
					alogpfrag[i] = 76;
					ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[76]);
					return;
				}
			}

			boolean flag1 = false;
			boolean flag2 = false;
			perturbacion = 0.0;
			if (fragment[i].equals("dNH"))
				perturbacion = refracval[50];
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					continue;
				if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
					if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
						alogpfrag[i] = 74;
						ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[74] + perturbacion);
						return;
					} else {
						flag1 = true;
					}
				} else if (!((IAtom)ca.get(j)).getSymbol().equals("C") || ((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
						flag2 = true;
				  else  flag2 = true;  // quitar esto *****Importante*******

				if (flag1 && flag2) { // X-N=X or Ar-N=X
					alogpfrag[i] = 78;
					ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[78] + perturbacion);
				} else {
					//logger.debug("missing group: R-N=X");
				}
			}
		} else if (fragment[i].indexOf("p") > -1)
			ai.setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[79] + perturbacion);
	}

	public void calcGroup081_to_085(int i) {

		if (!fragment[i].equals("sF"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
		int doublebondcount = 0;
		int triplebondcount = 0;
		String hybrid = "";
		for (int j = 0; j <= bonds.size() - 1; j++) {
			IBond bj = (IBond)bonds.get(j);
			if (bj.getOrder() == IBond.Order.DOUBLE)
				doublebondcount++;
			else if (bj.getOrder() == IBond.Order.TRIPLE)
				triplebondcount++;
		}
		if (doublebondcount == 0 && triplebondcount == 0)
			hybrid = "sp3";
		else if (doublebondcount == 1)
			hybrid = "sp2";
		else if (doublebondcount == 2 || triplebondcount == 1)
			hybrid = "sp";

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		int OxNum = 0;
		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1)
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
		}
		if (hybrid.equals("sp3") && OxNum == 1) {
			alogpfrag[i] = 81;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[81]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			alogpfrag[i] = 82;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[82]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			alogpfrag[i] = 83;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[83]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			alogpfrag[i] = 84;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[84]);
		} else if ((hybrid.equals("sp2") && OxNum > 1) || (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4) || !ca0.getSymbol().equals("C")) {
			alogpfrag[i] = 85;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[85]);
		}

	}

	public void calcGroup086_to_090(int i) {

		if (!fragment[i].equals("sCl"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
		int doublebondcount = 0;
		int triplebondcount = 0;
		String hybrid = "";
		for (int j = 0; j <= bonds.size() - 1; j++) {
			IBond bj = (IBond)bonds.get(j);
			if (bj.getOrder() == IBond.Order.DOUBLE)
				doublebondcount++;
			else if (bj.getOrder() == IBond.Order.TRIPLE)
				triplebondcount++;
		}
		if (doublebondcount == 0 && triplebondcount == 0)
			hybrid = "sp3";
		else if (doublebondcount == 1)
			hybrid = "sp2";
		else if (doublebondcount == 2 || triplebondcount == 1)
			hybrid = "sp";

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		int OxNum = 0;
		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			String s = ca2j.getSymbol();
			if (ap.getNormalizedElectronegativity(s) > 1)
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
		}
		if (hybrid.equals("sp3") && OxNum == 1) {
			alogpfrag[i] = 86;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[86]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			alogpfrag[i] = 87;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[87]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			alogpfrag[i] = 88;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[88]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			alogpfrag[i] = 89;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[89]);
		} else if ((hybrid.equals("sp2") && OxNum > 1) || (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4) || !ca0.getSymbol().equals("C")) {
			alogpfrag[i] = 90;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[90]);
		}
	}

	public void calcGroup091_to_095(int i) {

		if (!fragment[i].equals("sBr"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
		int doublebondcount = 0;
		int triplebondcount = 0;
		String hybrid = "";
		for (int j = 0; j <= bonds.size() - 1; j++) {
			IBond bj = (IBond)bonds.get(j);
			if (bj.getOrder() == IBond.Order.DOUBLE)
				doublebondcount++;
			else if (bj.getOrder() == IBond.Order.TRIPLE)
				triplebondcount++;
		}
		if (doublebondcount == 0 && triplebondcount == 0)
			hybrid = "sp3";
		else if (doublebondcount == 1)
			hybrid = "sp2";
		else if (doublebondcount == 2 || triplebondcount == 1)
			hybrid = "sp";

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		int OxNum = 0;
		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1)
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
		}

		if (hybrid.equals("sp3") && OxNum == 1) {
			alogpfrag[i] = 91;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[91]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			alogpfrag[i] = 92;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[92]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			alogpfrag[i] = 93;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[93]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			alogpfrag[i] = 94;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[94]);
		} else if ((hybrid.equals("sp2") && OxNum > 1) || (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4) || !ca0.getSymbol().equals("C")) {
			alogpfrag[i] = 95;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[95]);
		}
	}

	public void calcGroup096_to_100(int i) {

		if (!fragment[i].equals("sI"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
		int doublebondcount = 0;
		int triplebondcount = 0;
		String hybrid = "";
		for (int j = 0; j <= bonds.size() - 1; j++) {
			IBond bj = (IBond)bonds.get(j);
			if (bj.getOrder() == IBond.Order.DOUBLE)
				doublebondcount++;
			else if (bj.getOrder() == IBond.Order.TRIPLE)
				triplebondcount++;
		}

		if (doublebondcount == 0 && triplebondcount == 0)
			hybrid = "sp3";
		else if (doublebondcount == 1)
			hybrid = "sp2";
		else if (doublebondcount == 2 || triplebondcount == 1)
			hybrid = "sp";

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		int OxNum = 0;
		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1)
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
		}

		if (hybrid.equals("sp3") && OxNum == 1) {
			alogpfrag[i] = 96;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[96]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
			alogpfrag[i] = 97;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[97]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
			alogpfrag[i] = 98;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[98]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
			alogpfrag[i] = 99;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[99]);
		} else if ((hybrid.equals("sp2") && OxNum > 1) || (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4) || !ca0.getSymbol().equals("C")) {
			alogpfrag[i] = 100;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[100]);
		}
	}

	public void calcGroup101_to_104(int i) {
		IAtom ai = atomContainer.getAtom(i);
		String s = ai.getSymbol();

		if (ai.getFormalCharge() == -1) {
			if (s.equals("F")) {
				alogpfrag[i] = 101;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[101]);
			} else if (s.equals("Cl")) {
				alogpfrag[i] = 102;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[102]);
			} else if (s.equals("Br")) {
				alogpfrag[i] = 103;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[103]);
			} else if (s.equals("I")) {
				alogpfrag[i] = 104;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[104]);
			}
		}
	}

	public void calcGroup106(int i) {
		if (fragment[i].equals("sSH")) {
			alogpfrag[i] = 106;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[106] + refracval[50]);
		}
	}

	public void calcGroup107(int i) {
		// S in R2S, RS-SR
		// R = any group linked through C
		// if (!Fragment[i].equals("SssS")) return;

		// In ALOGP, for malathion PSC is consider to have group 107 (even
		// though has P instead of R)

		// for lack of fragment, use this fragment for SaaS
        
		//Quitar sSm *******Importante*********
		if (fragment[i].equals("ssS") || fragment[i].equals("aaS") || fragment[i].equals("sSm")) {
			alogpfrag[i] = 107;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[107]);
		}
		// IAtom [] ca=atomContainer.getConnectedAtoms(atomContainer.getAtomAt(i));
		//
		// if ((ca[0].getSymbol().equals("C") && ca[1].getSymbol().equals("C"))
		// ||
		// (ca[0].getSymbol().equals("C") && ca[1].getSymbol().equals("S")) ||
		// (ca[0].getSymbol().equals("S") && ca[1].getSymbol().equals("C"))) {
		// frags[107]++;
		// alogpfrag[i]=107;
		// }
	}

	public void calcGroup108(int i) {
		if (fragment[i].equals("dS")) {
			alogpfrag[i] = 108;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[108]);
		}
	}

	public void calcGroup109(int i) {
		if (!fragment[i].equals("dssS")) return;
		
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);
		int SdOCount=0;
		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
				}
			} else if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("O"))
					SdOCount++;
			}
		}
		if (SdOCount==1) { // for now dont check if SsCCount==2
			alogpfrag[i] = 109;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[109]);
		}
	}

	public void calcGroup110(int i) {
		if (!fragment[i].equals("ddssS"))
			return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);
		int SdOCount=0;
		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
				}
			} else if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
				//Quitar ((IAtom)ca.get(j)).getSymbol().equals("N") ***********Importante***********
				if (((IAtom)ca.get(j)).getSymbol().equals("O") || ((IAtom)ca.get(j)).getSymbol().equals("N"))
					SdOCount++;
			}
		}
		
		//Quitar SdOCount==1 *****importante*******
		if (SdOCount==2 || SdOCount==1) { // for now dont check if SsCCount==2
			alogpfrag[i] = 110;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[110]);
		}

	}

	public void calcGroup111(int i) {
		if (fragment[i].equals("ssssSi")) {
			alogpfrag[i] = 111;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[111]);
		}
	}

	public void calcGroup116_117_120(int i) {
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);
		int XCount=0;
		int RCount=0;
		boolean PdX=false;
		if (!fragment[i].equals("dsssP")) return;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					RCount++;
				else
					XCount++;
			} else if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C"))
					PdX=true;
			}
		}
		if (PdX) {
			if (RCount == 3) {
				alogpfrag[i] = 116;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[116]);
			} else if (XCount == 3) {
				alogpfrag[i] = 117;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[117]);
			} else if (XCount == 2 && RCount == 1) {
				alogpfrag[i] = 120;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[120]);
			}
			
			//quitar esto **********importante************
			else if (XCount == 1 && RCount == 2) {
				alogpfrag[i] = 120;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[120]);
			}
		}
	}
	
    public void calcGroup118_119(int i) {
		if (!fragment[i].equals("sssP")) return;

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);
		int XCount=0;
		int RCount=0;
		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					RCount++;
				else
					XCount++;
			}
		}
		if (XCount==3) {
			alogpfrag[i] = 118;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[118]);
		} else if (RCount==3) {
			alogpfrag[i] = 119;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE, refracval[119]);
		}
	}

	public DescriptorValue calculate(IAtomContainer atomContainer) {
		try {
			container = (IAtomContainer) atomContainer.clone();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);            
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(container.getBuilder());
			hAdder.addImplicitHydrogens(container);
			//CDKHueckelAromaticityDetector.detectAromaticity(container);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
			boolean aromatic = CDKHueckelAromaticityDetector.detectAromaticity(container);
			container.setFlag(CDKConstants.ISAROMATIC, aromatic);
			this.atomContainer = container;
		} catch (CloneNotSupportedException e) {
			return getDummyDescriptorValue(new CDKException("Error during clone"));
		} catch (CDKException e) {
			return getDummyDescriptorValue(new CDKException("Error during atom typing" + e.getMessage()));
		}

		IRingSet rs;
		try {
			AllRingsFinder arf = new AllRingsFinder();
			rs = arf.findAllRings(container);
			//Mio
			selectAromaticRing(rs);

		} catch (Exception e) {
			return getDummyDescriptorValue(new CDKException("Could not find all rings: " + e.getMessage()));
		}


		String[] fragment = new String[container.getAtomCount()];
		EStateAtomTypeMatcher eStateMatcher = new EStateAtomTypeMatcher();
		eStateMatcher.setRingSet(rs);

		for (int i = 0; i < container.getAtomCount(); i++) {
			IAtomType atomType = eStateMatcher.findMatchingAtomType(container, container.getAtom(i));			
			if (atomType == null)
				fragment[i] = null;
			else 
				fragment[i] = atomType.getAtomTypeName();
		}

		try {
		  if(!isRIntrinsicValue())
			calculate(container, fragment, rs);
		} catch (CDKException e) {
			return getDummyDescriptorValue(new CDKException(e.getMessage()));
		}

		double total = Double.NaN;
		try {
			total = updateRefractedTopographicValues();
			this.container.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC, total);
		} catch (CDKException e) {
			
			return getDummyDescriptorValue(new CDKException(e.getMessage()));
		}
		return new DescriptorValue(getSpecification(), getParameterNames(),
				getParameters(), new DoubleResult(total), getDescriptorNames());
	}

	public double updateRefractedTopographicValues() throws CDKException {

		IAtomicDescriptor distance = new DistanceToAtomDescriptor();
		double efecto = 0.0, total = 0.0;

		int nAtomos = container.getAtomCount();
		for (int i = 0; i < nAtomos; i++) {
			efecto = 0.0;
			IAtom atomoOrigen = container.getAtom(i);
			if(!atomoOrigen.getSymbol().equals("H")){
				double valorEstateIntAtomoOrigen = Double.valueOf(atomoOrigen.getProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE).toString());
				for(int j = 0; j < nAtomos; j++){
					IAtom atomoDestino = container.getAtom(j);
					if(i != j && !atomoDestino.getSymbol().equals("H")){
						Object[] objs = {j};
						distance.setParameters(objs);
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
	
	public boolean isRIntrinsicValue(){	    	
	    	boolean aux = false;
	    	try{
	    		//Primer atomo, si esta bien formado no debe ser H
	    		this.container.getFirstAtom().getProperty(IHybridDescriptor.Descriptor.REFRACTIINTRINSICVALUE).toString();
	    		aux = true;
	    	}catch (Exception e) {	
			}
	    	return aux;
	    }

}
