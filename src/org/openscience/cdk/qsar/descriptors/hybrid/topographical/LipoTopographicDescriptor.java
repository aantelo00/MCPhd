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

public class LipoTopographicDescriptor  extends TopographicDescriptor {

	final static double[] fragval = new double[121];// coefficients for alogp model
	static {
//		// fragments for ALOGP from Ghose et al., 1998
		fragval[1] = -1.5603;
		fragval[2] = -1.012;
		fragval[3] = -0.6681;
		fragval[4] = -0.3698;
		fragval[5] = -1.788;
		fragval[6] = -1.2486;
		fragval[7] = -1.0305;
		fragval[8] = -0.6805;
		fragval[9] = -0.3858;
		fragval[10] = 0.7555;
		fragval[11] = -0.2849;
		fragval[12] = 0.02;
		fragval[13] = 0.7894;
		fragval[14] = 1.6422;
		fragval[15] = -0.7866;
		fragval[16] = -0.3962;
		fragval[17] = 0.0383;
		fragval[18] = -0.8051;
		fragval[19] = -0.2129;
		fragval[20] = 0.2432;
		fragval[21] = 0.4697;
		fragval[22] = 0.2952;
		fragval[23] = 0;
		fragval[24] = -0.3251;
		fragval[25] = 0.1492;
		fragval[26] = 0.1539;
		fragval[27] = 0.0005;
		fragval[28] = 0.2361;
		fragval[29] = 0.3514;
		fragval[30] = 0.1814;
		fragval[31] = 0.0901;
		fragval[32] = 0.5142;
		fragval[33] = -0.3723;
		fragval[34] = 0.2813;
		fragval[35] = 0.1191;
		fragval[36] = -0.132;
		fragval[37] = -0.0244;
		fragval[38] = -0.2405;
		fragval[39] = -0.0909;
		fragval[40] = -0.1002;
		fragval[41] = 0.4182;
		fragval[42] = -0.2147;
		fragval[43] = -0.0009;
		fragval[44] = 0.1388;
		fragval[45] = 0;
		fragval[46] = 0.7341;
		fragval[47] = 0.6301;
		fragval[48] = 0.518;
		fragval[49] = -0.0371;
		fragval[50] = -0.1036;
		fragval[51] = 0.5234;
		fragval[52] = 0.6666;
		fragval[53] = 0.5372;
		fragval[54] = 0.6338;
		fragval[55] = 0.362;
		fragval[56] = -0.3567;
		fragval[57] = -0.0127;
		fragval[58] = -0.0233;
		fragval[59] = -0.1541;
		fragval[60] = 0.0324;
		fragval[61] = 1.052;
		fragval[62] = -0.7941;
		fragval[63] = 0.4165;
		fragval[64] = 0.6601;
		fragval[65] = 0;
		fragval[66] = -0.5427;
		fragval[67] = -0.3168;
		fragval[68] = 0.0132;
		fragval[69] = -0.3883;
		fragval[70] = -0.0389;
		fragval[71] = 0.1087;
		fragval[72] = -0.5113;
		fragval[73] = 0.1259;
		fragval[74] = 0.1349;
		fragval[75] = -0.1624;
		fragval[76] = -2.0585;
		fragval[77] = -1.915;
		fragval[78] = 0.4208;
		fragval[79] = -1.4439;
		fragval[80] = 0;
		fragval[81] = 0.4797;
		fragval[82] = 0.2358;
		fragval[83] = 0.1029;
		fragval[84] = 0.3566;
		fragval[85] = 0.1988;
		fragval[86] = 0.7443;
		fragval[87] = 0.5337;
		fragval[88] = 0.2996;
		fragval[89] = 0.8155;
		fragval[90] = 0.4856;
		fragval[91] = 0.8888;
		fragval[92] = 0.7452;
		fragval[93] = 0.5034;
		fragval[94] = 0.8995;
		fragval[95] = 0.5946;
		fragval[96] = 1.4201;
		fragval[97] = 1.1472;
		fragval[98] = 0;
		fragval[99] = 0.7293;
		fragval[100] = 0.7173;
		fragval[101] = 0;
		fragval[102] = -2.6737;
		fragval[103] = -2.4178;
		fragval[104] = -3.1121;
		fragval[105] = 0;
		fragval[106] = 0.6146;
		fragval[107] = 0.5906;
		fragval[108] = 0.8758;
		fragval[109] = -0.4979;
		fragval[110] = -0.3786;
		fragval[111] = 1.5188;
		fragval[112] = 1.0255;
		fragval[113] = 0;
		fragval[114] = 0;
		fragval[115] = 0;
		fragval[116] = -0.9359;
		fragval[117] = -0.1726;
		fragval[118] = -0.7966;
		fragval[119] = 0.6705;
		fragval[120] = -0.4801;
	}
	
	public LipoTopographicDescriptor() throws CDKException {
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
		if (fragment[i].equals("sCH3")) {
			perturbacion = 0.0;		
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j < ca.size(); j++)
			{
				IAtom atomoAux = ca.get(j);
				if (atomoAux.getSymbol().equals("C")) {
					alogpfrag[i] = 1;
				}
				/*IAtom atomoAux = ca.get(j);
				if (atomoAux.equals("C")) {
					alogpfrag[i] = 1;
				}*/
				else if (atomoAux.getSymbol().equals("H")) {
					perturbacion += fragval[htype];
				}
				else {
					alogpfrag[i] = 5;
					atomos5++;
				}
		   }	
		   if(atomos5 == 0)
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[1] + perturbacion);
		   else 
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[5] + perturbacion);
		   return;
		}
		if(fragment[i].equals("ssCH2")) {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j < ca.size(); j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else if (((IAtom)ca.get(j)).getSymbol().equals("H")) 
				    perturbacion += fragval[htype];
				else
					HeteroCount++;
		    }
			if (CarbonCount == 2 && HeteroCount == 0) {
				alogpfrag[i] = 2;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[2] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 1) {
				alogpfrag[i] = 6;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[6] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 2) {
				alogpfrag[i] = 7;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[7] + perturbacion);
			}
			return;
		}
		if(fragment[i].equals("sssCH")) {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					perturbacion += fragval[htype];
			    else
					HeteroCount++;
			}
			if (CarbonCount == 3 && HeteroCount == 0) {
				alogpfrag[i] = 3;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[3] + perturbacion);
			} else if (CarbonCount == 2 && HeteroCount == 1) {
				alogpfrag[i] = 8;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[8] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 2) {
				alogpfrag[i] = 9;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[9] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 3) {
				alogpfrag[i] = 10;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[10] + perturbacion);
			}
			return;	
		}
		if(fragment[i].equals("ssssC")) {
			perturbacion = 0.0;
			CarbonCount = 0;
			HeteroCount = 0;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C"))
					CarbonCount++;
				else
					HeteroCount++;
	        }
			if (CarbonCount == 4 && HeteroCount == 0) {
				alogpfrag[i] = 4;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[4] + perturbacion);
			} else if (CarbonCount == 3 && HeteroCount == 1) {
				alogpfrag[i] = 11;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[11] + perturbacion);
			} else if (CarbonCount == 2 && HeteroCount == 2) {
				alogpfrag[i] = 12;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[12] + perturbacion);
			} else if (CarbonCount == 1 && HeteroCount == 3) {
				alogpfrag[i] = 13;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[13] + perturbacion);
			} else if (CarbonCount == 0 && HeteroCount == 4) {
				alogpfrag[i] = 14;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[14] + perturbacion);
			}
			return;	
		}
		if(fragment[i].equals("dCH2")) {
			perturbacion = 0.0;
			htype = getHAtomType(atomContainer.getAtom(i), null);
			alogpfrag[i] = 15;
			perturbacion += 2 * fragval[htype];
			atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[15] + perturbacion);
            return;		
		}
		if(fragment[i].equals("dsCH")) {
			perturbacion = 0.0;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			perturbacion += fragval[htype];
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (((IAtom)ca.get(j)).getSymbol().equals("H")){
					continue;
				}
				if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C")) 
						HaveCsX = true;
					if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
						HaveCsAr = true;
				} else if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C"))
						HaveCdX = true;
				}
			}
			if (HaveCdX) {
				if (HaveCsAr) {
					alogpfrag[i] = 37;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[37] + perturbacion);
				} else {
					alogpfrag[i] = 36;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[36] + perturbacion);
				}
			} else {
				if (HaveCsX) {
					alogpfrag[i] = 18;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[18] + perturbacion);
				} else {
					alogpfrag[i] = 16;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[16] + perturbacion);
				}
			}
			return;
		}
		if(fragment[i].equals("dssC")) {
			HaveCdX = false;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
					if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					} else {
					}
					if (!((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC)) {
					} else {
					}
				} else if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
					if (!((IAtom)ca.get(j)).getSymbol().equals("C"))
						HaveCdX = true;
				}
			}
			if (HaveCdX) {
				if (HaveCsAr) {
					alogpfrag[i] = 37;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[37] + perturbacion);
				} else {
					alogpfrag[i] = 36;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[36] + perturbacion);
				}
			} else {
				if (HaveCsX) {
					alogpfrag[i] = 18;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[18] + perturbacion);
				} else {
					alogpfrag[i] = 16;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[16] + perturbacion);
				}
			}
			return;	
		}
		if(fragment[i].equals("tCH")) {
			alogpfrag[i] = 21;
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[21] + fragval[htype]);
		    return;
		}
		if(fragment[i].equals("ddC")) {
			if (((IAtom)ca.get(0)).getSymbol().equals("C") && ((IAtom)ca.get(1)).getSymbol().equals("C")){
				alogpfrag[i] = 22;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[22]);
			} else if (!((IAtom)ca.get(0)).getSymbol().equals("C") && !((IAtom)ca.get(1)).getSymbol().equals("C")) {
				alogpfrag[i] = 40;
			    atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[40]);
			}
			return;
		}
		if(fragment[i].equals("tsC")) {
			HaveCtX = false;
			HaveCsX = false;
			for (int j = 0; j <= ca.size() - 1; j++) {
				if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE)
					if (!((IAtom)ca.get(j)).getSymbol().equals("C"))
						HaveCsX = true;
				else if (atomContainer.getBond(atomoPrincipal, ((IAtom)ca.get(j))).getOrder() == IBond.Order.TRIPLE)
					if (!((IAtom)ca.get(j)).getSymbol().equals("C"))
						HaveCtX = true;
			}
			if (HaveCtX && !HaveCsX) {
				alogpfrag[i] = 40;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[40]);
			} else if (HaveCsX) {// #C-X
				alogpfrag[i] = 23;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[23]);
			} else if (!HaveCsX) { // #C-R
				alogpfrag[i] = 22;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[22]);
			}
			return;
		}
		if (fragment[i].equals("aaCH")) {
			htype = getHAtomType(atomContainer.getAtom(i), ca);
			IAtom ca0;
			IAtom ca1;
			if (((IAtom)ca.get(0)).getSymbol().equals("H")) {
				ca0 = (IAtom)ca.get(1);
				ca1 = (IAtom)ca.get(2);
			} else {
				if (((IAtom)ca.get(1)).getSymbol().equals("H")) {
					ca0 = (IAtom)ca.get(0);
					ca1 = (IAtom)ca.get(2);
				} else {
					ca0 = (IAtom)ca.get(0);
					ca1 = (IAtom)ca.get(1);
				}
			}
			if (ca0.getSymbol().equals("C") && ca1.getSymbol().equals("C")) {
				alogpfrag[i] = 24;
				atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[24] + fragval[htype]);
				return;
			}
			List<IBond> bonds = atomContainer.getConnectedBondsList(ca0);
			boolean HaveDouble1 = false;
			for (int k = 0; k <= bonds.size() - 1; k++)
			{
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble1 = true;
					break;
				}
			}
			bonds = atomContainer.getConnectedBondsList(ca1);
			boolean HaveDouble2 = false;
			for (int k = 0; k <= bonds.size() - 1; k++) {
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble2 = true;
					break;
				}
			}
			if (!(ca0).getSymbol().equals("C") && !((IAtom)ca.get(1)).getSymbol().equals("C")) {
				if (HaveDouble1 && HaveDouble2) { // X--CH--X
					alogpfrag[i] = 30;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[30] + fragval[htype]);
				} else { // X--CH...X
					alogpfrag[i] = 42;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[42] + fragval[htype]);
				}

			} else if (ca0.getSymbol().equals("C") && !ca1.getSymbol().equals("C")
					|| (!ca0.getSymbol().equals("C") && ca1.getSymbol().equals("C"))) {
				if (HaveDouble1 && HaveDouble2) { // R--CH--X
					alogpfrag[i] = 27;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[27] + fragval[htype]);
				} else {// R--CH...X
					alogpfrag[i] = 33;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[33] + fragval[htype]);
				}
			}
            return;
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
			for (int k = 0; k <= bonds.size() - 1; k++) {
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble1 = true;
					break;
				}
			}
			bonds = atomContainer.getConnectedBondsList(sameringatoms[1]);
			boolean HaveDouble2 = false;
			for (int k = 0; k <= bonds.size() - 1; k++) {
				if (((IBond)bonds.get(k)).getOrder() == IBond.Order.DOUBLE) {
					HaveDouble2 = true;
					break;
				}
			}
			if (!sameringatoms[0].getSymbol().equals("C")
					&& !sameringatoms[1].getSymbol().equals("C")) {
				if (HaveDouble1 && HaveDouble2) { // X--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						alogpfrag[i] = 31;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[31]);
					} else { // X--CX--X
						alogpfrag[i] = 32;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[32]);
					}
				} else {
					if (nonringatom.getSymbol().equals("C")) { // X--CR..X
						alogpfrag[i] = 43;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[43]);
					} else { // X--CX...X
						alogpfrag[i] = 44;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[44]);
					}
				}
			} else if (sameringatoms[0].getSymbol().equals("C")
					&& sameringatoms[1].getSymbol().equals("C")) {
				if (nonringatom.getSymbol().equals("C")) {// R--CR--R
					alogpfrag[i] = 25;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[25]);
				} else { // R--CX--R
					alogpfrag[i] = 26;
					atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[26]);
				}
			} else if ((sameringatoms[0].getSymbol().equals("C") && !sameringatoms[1].getSymbol().equals("C"))
			           || (!sameringatoms[0].getSymbol().equals("C") && sameringatoms[1].getSymbol().equals("C"))) {
				if (HaveDouble1 && HaveDouble2) { // R--CR--X
					if (nonringatom.getSymbol().equals("C")) {
						alogpfrag[i] = 28;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[28]);
					} else { // R--CX--X
						alogpfrag[i] = 29;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[29]);
					}
				} else {
					if (nonringatom.getSymbol().equals("C")) { // R--CR..X
						alogpfrag[i] = 34;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[34]);
					} else { // R--CX...X
						alogpfrag[i] = 35;
						atomoPrincipal.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[35]);
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
//		frags[50]++; //H atom attached to a hetero atom

		IAtom ca0 = (IAtom)ca.get(0);
		if (ca0.getSymbol().equals("H"))
			ca0 = (IAtom)ca.get(1);

		if (ca0.getFlag(CDKConstants.ISAROMATIC)) { // phenol
//			frags[57]++;
			alogpfrag[i] = 57;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[57] + fragval[50]);

			return;
		}

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);
		for (int j = 0; j <= ca2.size() - 1; j++) {
			if (atomContainer.getBond((IAtom)ca2.get(j), ca0).getOrder() == IBond.Order.DOUBLE) {
//				frags[57]++;
				alogpfrag[i] = 57;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[57] + fragval[50]);
				return;
			}
		}
//		frags[56]++;
		alogpfrag[i] = 56;
		atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[56] + fragval[50]);
	}

	public void calcGroup058_61(int i) {
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));

		// 58: O in =O
		// 61: --O in nitro, N-oxides
		// 62: O in O-
		IAtom ca0 = (IAtom)ca.get(0);

		if (fragment[i].equals("sOm")) {

			if (ca0.getSymbol().equals("N") && ca0.getFormalCharge() == 1) {
//				frags[61]++;
				alogpfrag[i] = 61;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[61]);
			} else {
//				frags[62]++;
				alogpfrag[i] = 62;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[62]);
			}

		} else if (fragment[i].equals("dO")) {
			if (ca0.getSymbol().equals("N") && ca0.getFormalCharge() == 1) {
//				frags[61]++;
				alogpfrag[i] = 61;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[61]);
			} else {
//				frags[58]++;
				alogpfrag[i] = 58;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[58]);
			}
		}

	}

	public void calcGroup059_060_063(int i) {
		// O in Al-O-Ar, Ar2O, R...O...R, ROC=X
		// ... = aromatic single bonds
		if (!fragment[i].equals("ssO") && !fragment[i].equals("aaO") && !fragment[i].equals("aaOp"))
			return;

		// Al-O-Ar, Ar2O
		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ca0 = (IAtom)ca.get(0);
		IAtom ca1 = (IAtom)ca.get(1);

		if (fragment[i].equals("ssO")) {
			if (ca0.getFlag(CDKConstants.ISAROMATIC)
					|| ca1.getFlag(CDKConstants.ISAROMATIC)) {
//				frags[60]++;
				alogpfrag[i] = 60;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[60]);

			} else {

				for (int j = 0; j <= ca.size() - 1; j++) {
					// if (((IAtom)ca.get(j)).getSymbol().equals("C")) { // for malathion
					// O-P(=S)
					// was considered to count as group 60

					List<IAtom> ca2 = atomContainer.getConnectedAtomsList(((IAtom)ca.get(j)));
					for (int k = 0; k <= ca2.size() - 1; k++) {
						if (atomContainer.getBond(((IAtom)ca.get(j)), (IAtom)ca2.get(k)).getOrder() == IBond.Order.DOUBLE) {
							if (!((IAtom)ca2.get(k)).getSymbol().equals("C")) {
//								frags[60]++;
//								alogpfrag[i] = 60;
								atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[60]);
								return;
							}
						}
					}

				} // end j ca loop

				if (ca0.getSymbol().equals("O")
						|| ca1.getSymbol().equals("O")) {
//					frags[63]++;
					alogpfrag[i] = 63;
					atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[63]);
				} else {
//					frags[59]++;
					alogpfrag[i] = 59;
					atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[59]);

				}

			}
		//Quitar aaOp  **********Importante************	
		} else if (fragment[i].equals("aaO") || fragment[i].equals("aaOp")) {
//			frags[60]++;
			alogpfrag[i] = 60;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[60]);
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
		//IAtom ca0 = (IAtom)ca.get(0);
		//IAtom ca1 = (IAtom)ca.get(1);

		for (int j = 0; j <= ca.size() - 1; j++)
		{
			if (((IAtom)ca.get(j)).getSymbol().equals("H"))
				continue;
			if (((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
				NAr++;
			else
				NAl++;
		}

		// first check if have RC(=O)N or NX=X
		for (int j = 0; j <= ca.size() - 1; j++)
		{
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
//								frags[72]++;
								alogpfrag[i] = 72;
								ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[72]);
								return;
							}
						}
					}
				}
			}
		}

		if (fragment[i].equals("sNH2"))
		{
			IAtom ca0 = null;
			//Find which neigbpur is not the hydrogen atom
			for (int j = 0; j <= ca.size() - 1; j++)
			{
				if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					continue;
				else
				{
					ca0 = (IAtom)ca.get(j);
					break;
				}
			}
			if (ca0.getFlag(CDKConstants.ISAROMATIC)
					|| !ca0.getSymbol().equals("C"))
			{
//				frags[69]++;
				alogpfrag[i] = 69;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[69] + 2 * fragval[50]);
			}
			else
			{
//				frags[66]++;
				alogpfrag[i] = 66;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[66] + 2 * fragval[50]);
			}
//			frags[50]+=2; //H atom attached to a hetero atom
		}
		else if (fragment[i].equals("aaNH") || fragment[i].equals("saaN"))
		{ // R...NH...R
			perturbacion = 0.0;
//			frags[73]++;
			alogpfrag[i] = 73;

			if (fragment[i].equals("aaNH")){
//				frags[50]++; //H atom attached to a hetero atom
				perturbacion = fragval[50];
			}
			ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[73] + perturbacion);
		}
		else if (fragment[i].equals("ssNH"))
		{
			if (NAr == 2 && NAl == 0) { // Ar2NH
//				frags[73]++;
				alogpfrag[i] = 73;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[73] + fragval[50]);

			} else if (NAr == 1 && NAl == 1) { // Ar-NH-Al
//				frags[70]++;
				alogpfrag[i] = 70;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[70] + fragval[50]);
			} else if (NAr == 0 && NAl == 2) { // Al2NH
//				frags[67]++;
				alogpfrag[i] = 67;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[67] + fragval[50]);
			}
//			frags[50]++; //H atom attached to a hetero atom
		}
		else if (fragment[i].equals("sssN"))
		{
			if ((NAr == 3 && NAl == 0) || (NAr == 2 && NAl == 1)) { // Ar3N &
				// Ar2NAl
//				frags[73]++;
				alogpfrag[i] = 73;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[73]);
			} else if (NAr == 1 && NAl == 2) {
//				frags[71]++;
				alogpfrag[i] = 71;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[71]);
			} else if (NAr == 0 && NAl == 3) {
//				frags[68]++;
				alogpfrag[i] = 68;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[68]);
			}
		}
		//arreglar esto, quitar aaaN y aaNm
		else if (fragment[i].equals("aaN") || fragment[i].equals("aaaN") || fragment[i].equals("aaNm") )
		{
//			frags[75]++;
			alogpfrag[i] = 75;
			ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[75]);
		}
		else if (fragment[i].equals("ssdNp") || fragment[i].equals("dssNp"))
		{
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
//				frags[76]++;
				alogpfrag[i] = 76;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[76]);
			} else if (HaveSsOm && HaveSdO && !Ar) {
//				frags[77]++;
				alogpfrag[i] = 77;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[77]);
			} else {
//				frags[79]++;
				alogpfrag[i] = 79;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[79]);
			}

		}
		else if (fragment[i].equals("tN"))
		{
			IAtom ca0 = (IAtom)ca.get(0);
			//Quitar  !ca0.getSymbol().equals("C") ***********Importante**************
			if (ca0.getSymbol().equals("C") || !ca0.getSymbol().equals("C")) { // R#N 
				alogpfrag[i] = 74;
				ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[74]);
			}
		}
		else if (fragment[i].equals("dNH") || fragment[i].equals("dsN"))
		{
			// test for RO-NO
			if (fragment[i].equals("dsN"))
			{
				IAtom ca0 = (IAtom)ca.get(0);
				IAtom ca1 = (IAtom)ca.get(1);
				if (ca0.getSymbol().equals("O")
						&& ca1.getSymbol().equals("O")) {
//					frags[76]++;
					alogpfrag[i] = 76;
					ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[76]);
					return;
				}
			}

			boolean flag1 = false;
			boolean flag2 = false;
			perturbacion = 0.0;
			if (fragment[i].equals("dNH")){
//				frags[50]++; //H atom attached to a hetero atom
				perturbacion = fragval[50];
			}
			for (int j = 0; j <= ca.size() - 1; j++)
			{
				if (((IAtom)ca.get(j)).getSymbol().equals("H"))
					continue;
				if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE)
				{
					if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
//						frags[74]++;
						alogpfrag[i] = 74;
						ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[74] + perturbacion);
						return;
					} else
						flag1 = true;
				} else if (!((IAtom)ca.get(j)).getSymbol().equals("C") || ((IAtom)ca.get(j)).getFlag(CDKConstants.ISAROMATIC))
						flag2 = true;
				  else  flag2 = true;  // quitar esto *****Importante*******

				if (flag1 && flag2)
				{ // X-N=X or Ar-N=X
//					frags[78]++;
					alogpfrag[i] = 78;
					ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[78] + perturbacion);
				} else
				{
					//logger.debug("missing group: R-N=X");
				}
			}

			//            if (fragment[i].equals("SdNH"))
			//                frags[50]++; //H atom attached to a hetero atom
		}
		else if (fragment[i].indexOf("p") > -1)
		{
//			frags[79]++;
//			alogpfrag[i] = 79;
			ai.setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[79] + perturbacion);
		}

		// TODO add code for R--N(--R)--O
		// first need to have program correctly read in structures with this
		// fragment (pyridine-n-oxides)
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
			if (bj.getOrder() == IBond.Order.DOUBLE) {
				doublebondcount++;
			}

			else if (bj.getOrder() == IBond.Order.TRIPLE) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);

		int OxNum = 0;

		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1) {
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
			}

		}

		if (hybrid.equals("sp3") && OxNum == 1) {
//			frags[81]++;
			alogpfrag[i] = 81;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[81]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
//			frags[82]++;
			alogpfrag[i] = 82;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[82]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
//			frags[83]++;
			alogpfrag[i] = 83;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[83]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
//			frags[84]++;
			alogpfrag[i] = 84;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[84]);
		} else if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !ca0.getSymbol().equals("C")) {
//			frags[85]++;
			alogpfrag[i] = 85;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[85]);
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
			if (bj.getOrder() == IBond.Order.DOUBLE) {
				doublebondcount++;
			}

			else if (bj.getOrder() == IBond.Order.TRIPLE) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);

		int OxNum = 0;

		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			String s = ca2j.getSymbol();

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (ap.getNormalizedElectronegativity(s) > 1) {
				// // F,O,Cl,Br,N
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
			}
		}

		if (hybrid.equals("sp3") && OxNum == 1) {
//			frags[86]++;
			alogpfrag[i] = 86;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[86]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
//			frags[87]++;
			alogpfrag[i] = 87;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[87]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
//			frags[88]++;
			alogpfrag[i] = 88;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[88]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
//			frags[89]++;
			alogpfrag[i] = 89;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[89]);
		} else if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !ca0.getSymbol().equals("C")) {
//			frags[90]++;
			alogpfrag[i] = 90;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[90]);
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
			if (bj.getOrder() == IBond.Order.DOUBLE) {
				doublebondcount++;
			}

			if (bj.getOrder() == IBond.Order.TRIPLE) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);

		int OxNum = 0;

		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1) {
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
			}

		}

		if (hybrid.equals("sp3") && OxNum == 1) {
//			frags[91]++;
			alogpfrag[i] = 91;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[91]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
//			frags[92]++;
			alogpfrag[i] = 92;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[92]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
//			frags[93]++;
			alogpfrag[i] = 93;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[93]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
//			frags[94]++;
			alogpfrag[i] = 94;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[94]);
		} else if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !ca0.getSymbol().equals("C")) {
//			frags[95]++;
			alogpfrag[i] = 95;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[95]);
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
			if (bj.getOrder() == IBond.Order.DOUBLE) {
				doublebondcount++;
			}

			else if (bj.getOrder() == IBond.Order.TRIPLE) {
				triplebondcount++;
			}

		}

		if (doublebondcount == 0 && triplebondcount == 0) {
			hybrid = "sp3";
		} else if (doublebondcount == 1) {
			hybrid = "sp2";
		} else if (doublebondcount == 2 || triplebondcount == 1) {
			hybrid = "sp";
		}

		List<IAtom> ca2 = atomContainer.getConnectedAtomsList(ca0);

		int OxNum = 0;

		for (int j = 0; j <= ca2.size() - 1; j++) {
			IAtom ca2j = (IAtom)ca2.get(j);
			

			// // F,O,Cl,Br,N

			// if (s.equals("F") || s.equals("O") || s.equals("Cl")
			// || s.equals("Br") || s.equals("N") || s.equals("S"))

			if (ap.getNormalizedElectronegativity(ca2j.getSymbol()) > 1) {
				OxNum += BondManipulator.destroyBondOrder(atomContainer.getBond(ca0, ca2j).getOrder());
			}

		}

		if (hybrid.equals("sp3") && OxNum == 1) {
//			frags[96]++;
			alogpfrag[i] = 96;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[96]);
		} else if (hybrid.equals("sp3") && OxNum == 2) {
//			frags[97]++;
			alogpfrag[i] = 97;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[97]);
		} else if (hybrid.equals("sp3") && OxNum == 3) {
//			frags[98]++;
			alogpfrag[i] = 98;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[98]);
		} else if (hybrid.equals("sp2") && OxNum == 1) {
//			frags[99]++;
			alogpfrag[i] = 99;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[99]);
		} else if ((hybrid.equals("sp2") && OxNum > 1)
				|| (hybrid.equals("sp") && OxNum >= 1)
				|| (hybrid.equals("sp3") && OxNum == 4)
				|| !ca0.getSymbol().equals("C")) {
//			frags[100]++;
			alogpfrag[i] = 100;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[100]);
		}

	}

	public void calcGroup101_to_104(int i) {
		IAtom ai = atomContainer.getAtom(i);
		String s = ai.getSymbol();

		if (ai.getFormalCharge() == -1) {
			if (s.equals("F")) {
//				frags[101]++;
				alogpfrag[i] = 101;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[101]);
			} else if (s.equals("Cl")) {
//				frags[102]++;
				alogpfrag[i] = 102;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[102]);
			} else if (s.equals("Br")) {
//				frags[103]++;
				alogpfrag[i] = 103;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[103]);
			} else if (s.equals("I")) {
//				frags[104]++;
				alogpfrag[i] = 104;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[104]);
			}

		}

	}

	public void calcGroup106(int i)
	{
		// S in SH
		if (fragment[i].equals("sSH")) {
//			frags[106]++;
			alogpfrag[i] = 106;
//			frags[50]++; //H atom attached to a hetero atom
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[106] + fragval[50]);
		}
	}

	public void calcGroup107(int i)
	{
		// S in R2S, RS-SR
		// R = any group linked through C
		// if (!Fragment[i].equals("SssS")) return;

		// In ALOGP, for malathion PSC is consider to have group 107 (even
		// though has P instead of R)

		// for lack of fragment, use this fragment for SaaS

		if (fragment[i].equals("ssS") || fragment[i].equals("aaS") || fragment[i].equals("sSm")) {
//			frags[107]++;
			alogpfrag[i] = 107;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[107]);
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

	public void calcGroup108(int i)
	{
		// S in R=S
		// In ALOGP, for malathion P=S is consider to have group 108 (even
		// though has P instead of R)
		if (fragment[i].equals("dS")) {
//			frags[108]++;
			alogpfrag[i] = 108;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[108]);
		}
	}

	public void calcGroup109(int i)
	{
		// for now S in O-S(=O)-O is assigned to this group
		// (it doesn't check which atoms are singly bonded to S
		if (!fragment[i].equals("dssS")) return;


		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);
		int SdOCount=0;
		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
				}
			} else if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("O")) {
					SdOCount++;
				}
			}
		}
		if (SdOCount==1) { // for now dont check if SsCCount==2
//			frags[109]++;
			alogpfrag[i] = 109;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[109]);
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
				if (((IAtom)ca.get(j)).getSymbol().equals("O") || ((IAtom)ca.get(j)).getSymbol().equals("N")) {
					SdOCount++;
				}
			}
		}
		if (SdOCount==2 || SdOCount==1) { // for now dont check if SsCCount==2
//			frags[110]++;
			alogpfrag[i] = 110;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[110]);
		}

	}

	public void calcGroup111(int i) {
		if (fragment[i].equals("ssssSi")) {
//			frags[111]++;
			alogpfrag[i] = 111;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[111]);
		}
	}

	public void calcGroup116_117_120(int i) {

		// S in R=S

		List<IAtom> ca = atomContainer.getConnectedAtomsList(atomContainer.getAtom(i));
		IAtom ai = atomContainer.getAtom(i);

		int XCount=0;
		int RCount=0;
		boolean PdX=false;

		if (!fragment[i].equals("dsssP")) return;

		for (int j = 0; j <= ca.size() - 1; j++) {
			if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.SINGLE) {
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					RCount++;
				} else {
					XCount++;
				}
			} else if (atomContainer.getBond(ai, ((IAtom)ca.get(j))).getOrder() == IBond.Order.DOUBLE) {
				if (!((IAtom)ca.get(j)).getSymbol().equals("C")) {
					PdX=true;
				}
			}
		}

		if (PdX) {
			if (RCount == 3) {
				alogpfrag[i] = 116;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[116]);
			} else if (XCount == 3) {
				alogpfrag[i] = 117;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[117]);
			} else if (XCount == 2 && RCount == 1) {
				alogpfrag[i] = 120;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[120]);
			}
			//Quitar esto **********importante*************
			else if (XCount == 1 && RCount == 2) {
				alogpfrag[i] = 120;
				atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[120]);
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
				if (((IAtom)ca.get(j)).getSymbol().equals("C")) {
					RCount++;
				} else {
					XCount++;
				}
			}
		}

		if (XCount==3) {
//			frags[118]++;
			alogpfrag[i] = 118;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[118]);
		} else if (RCount==3) {
//			frags[119]++;
			alogpfrag[i] = 119;
			atomContainer.getAtom(i).setProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE, fragval[119]);
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
		
//		Object value = molecule.getProperty("myProperty");
		


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
			//Comprobar el tipo de �tomo devuelto
			// System.out.println(atomType.toString());

			if (atomType == null) {
				fragment[i] = null;
			} else {
				fragment[i] = atomType.getAtomTypeName();
			}
		}

//		double[] ret = new double[0];
		try {
//			ret = calculate(container, fragment, rs);
			if(!isLIntrinsicValue())
			   calculate(container, fragment, rs);
		} catch (CDKException e) {
			return getDummyDescriptorValue(new CDKException(e.getMessage()));
		}

//		DoubleArrayResult results = new DoubleArrayResult();
//		results.add(ret[0]);
//		results.add(ret[1]);
//		results.add(ret[2]);

		double total = Double.NaN;
		try {
			total = updateLipoTopographicValues();
			this.container.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC, total);
		} catch (CDKException e) {
			
			return getDummyDescriptorValue(new CDKException(e.getMessage()));
		}

//		return new DescriptorValue(getSpecification(), getParameterNames(),
//				getParameters(), results, getDescriptorNames());
		return new DescriptorValue(getSpecification(), getParameterNames(),
				getParameters(), new DoubleResult(total), getDescriptorNames());
	}

	private double updateLipoTopographicValues() throws CDKException {

		IAtomicDescriptor distance = new DistanceToAtomDescriptor();
		double efecto = 0.0, total = 0.0;
		//Object[] arrAux = new Object[1];

		int nAtomos = container.getAtomCount();
		for (int i = 0; i < nAtomos; i++) {
			efecto = 0.0;
			IAtom atomoOrigen = container.getAtom(i);
			if(!atomoOrigen.getSymbol().equals("H")){
				double valorEstateIntAtomoOrigen = Double.valueOf(atomoOrigen.getProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE).toString());
				for(int j = 0; j < nAtomos; j++){
					IAtom atomoDestino = container.getAtom(j);
					if(i != j && !atomoDestino.getSymbol().equals("H")){
						Object[] objs = {j};
						distance.setParameters(objs);
						DescriptorValue value =  distance.calculate(atomoOrigen, container);
						double dist = ((DoubleResult) value.getValue()).doubleValue();
						if(dist != 0){
							double valorEstateIntAtomoDestino = Double.valueOf(atomoDestino.getProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE).toString());
							efecto += (valorEstateIntAtomoOrigen - valorEstateIntAtomoDestino)/Math.pow(dist, 2);
						}
					}
				}

				atomoOrigen.setProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC, valorEstateIntAtomoOrigen + efecto);
				total += valorEstateIntAtomoOrigen + efecto;

			}
		}
		return total;
	}
	
	private boolean isLIntrinsicValue(){

		boolean aux = false;
		try{
			//Primer atomo, si esta bien formado no debe ser H
			this.container.getFirstAtom().getProperty(IHybridDescriptor.Descriptor.LIPOINTRINSICVALUE).toString();
			aux = true;

		}catch (Exception e) {

		}

		return aux;
	}



}
