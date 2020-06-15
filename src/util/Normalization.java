package util;

import java.text.NumberFormat;
import java.util.ArrayList;

import org.openscience.cdk.qsar.IHybridDescriptor;

import domain.Molecule;

public class Normalization {

	private double vMin;
	private double vMaxList;
	private NumberFormat format;
	
	public Normalization(Molecule molecule) {
		format = NumberFormat.getInstance();
		format.setMaximumFractionDigits(0);
		ArrayList<Double> listValues = new ArrayList<Double>();
		listValues = newList(molecule);
		vMin = getMenor(listValues);
		ArrayList<Double> listValuesNormalization = new ArrayList<Double>();
		listValuesNormalization = normalizacionLista(molecule);
		@SuppressWarnings("unused")
		ArrayList<Double> listProportion = new ArrayList<Double>();
		listProportion = newListProportion(listValuesNormalization);
		vMaxList = getMayor(listValuesNormalization);
	}

	public ArrayList<Double> newList(Molecule molecule){
		ArrayList<Double> list = new ArrayList<Double>();
		for (int i=0; i< molecule.getCantAtmos(); i++) {
			if (!molecule.getGraph().getAtom(i).getSymbol().equals("H")) {
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC));
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC));
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC));
			}
		}
		return list;
	}
	
	public ArrayList<Double> newListProportion(ArrayList<Double> listValues){
		ArrayList<Double> list = new ArrayList<Double>();
		double mayor = getMayor(listValues);
		for(int i=0; i< listValues.size(); i ++)
		  list.add(listValues.get(i)/mayor);
		return list;
	}

	private ArrayList<Double> normalizacionLista(Molecule molecule) {
		ArrayList<Double> list = new ArrayList<Double>();
		for (int i = 0; i < molecule.getCantAtmos(); i++) {
			if (!molecule.getGraph().getAtom(i).getSymbol().equals("H")) {
                if(vMin > 0){
				    list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol()));
			        list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol()));
	                list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol()));
         	                 
		           }
                else {
				    list.add(((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + (-1 * vMin) + 1) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol()));
		            list.add(((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) + (-1 * vMin) + 1) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol()));
		            list.add(((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) + (-1 * vMin) + 1) / getRadioAtomico(molecule.getGraph().getAtom(i).getSymbol())); 
				}	
			}
		}
		return list;
	}

	public String getNormalization(double valor, String symbol) {
		double normal = normalizacionEstandar(valor, symbol, vMin);
		String rango = normalizarRango(normal, vMaxList);
		return rango;
	}

	public double normalizacionEstandar(double valor, String symbol, double vMinimo) {
		double ratio;
		if (vMinimo > 0)
			ratio = valor / getRadioAtomico(symbol);
		else
			ratio = (valor + (-1 * vMinimo) +1 ) / getRadioAtomico(symbol);
		return ratio;
	}

	public String normalizarRango(double valor, double vMaximo) {
		return format.format((valor/vMaximo)*100);
	}

	private double getMenor(ArrayList<Double> list) {
		double menor = list.get(0);
		for (int i = 1; i < list.size(); i++)
			if (menor > list.get(i))
				menor = list.get(i);
		return menor;
	}

	private double getMayor(ArrayList<Double> lista) {
		double mayor = lista.get(0);
		for (int i = 1; i < lista.size(); i++) 
			if (mayor < lista.get(i))
				mayor = lista.get(i);
		return mayor;
	}

	private double getRadioAtomico(String symbol) {
		Double valor = 0.0;
		if (symbol.equals("H"))
			valor = 0.23;
		else if (symbol.equals("He"))
			valor = 0.93;
		else if (symbol.equals("Li"))
			valor = 0.68;
		else if (symbol.equals("Be"))
			valor = 0.35;
		else if (symbol.equals("B"))
			valor = 0.83;
		else if (symbol.equals("C") || symbol.equals("N") || symbol.equals("O"))
			valor = 0.68;
		else if (symbol.equals("F"))
			valor = 0.64;
		else if (symbol.equals("Ne") || symbol.equals("Sr"))
			valor = 1.12;
		else if (symbol.equals("Na"))
			valor = 0.97;
		else if (symbol.equals("Mg"))
			valor = 1.1;
		else if (symbol.equals("Al") || symbol.equals("Cr") || symbol.equals("Mn") || symbol.equals("Tc") || symbol.equals("Re"))
			valor = 1.35;
		else if (symbol.equals("Si"))
			valor = 1.2;
		else if (symbol.equals("P"))
			valor = 0.75;
		else if (symbol.equals("S"))
			valor = 1.02;
		else if (symbol.equals("Cl") || symbol.equals("Ca"))
			valor = 0.99;
		else if (symbol.equals("Ar"))
			valor = 1.57;
		else if (symbol.equals("K") || symbol.equals("Co") || symbol.equals("V"))
			valor = 1.33;
		else if (symbol.equals("Sc"))
			valor = 1.44;
		else if (symbol.equals("Ti") || symbol.equals("Rb") || symbol.equals("Mo") || symbol.equals("Te"))
			valor = 1.47;
		else if (symbol.equals("Fe") || symbol.equals("Ba"))
			valor = 1.34;
		else if (symbol.equals("Ni") || symbol.equals("Pd") || symbol.equals("Pt") || symbol.equals("Au") || symbol.equals("Cm") || symbol.equals("Bk") || symbol.equals("Cf") || symbol.equals("Es") || symbol.equals("Fm") || symbol.equals("Md") || symbol.equals("No") || symbol.equals("Lr"))
			valor = 1.5;
		else if (symbol.equals("Cu"))
			valor = 1.52;
		else if (symbol.equals("Zn") || symbol.equals("Rh"))
			valor = 1.45;
		else if (symbol.equals("Ga") || symbol.equals("Se"))
			valor = 1.22;
		else if (symbol.equals("Ge"))
			valor = 1.17;
		else if (symbol.equals("As") || symbol.equals("Br"))
			valor = 1.21;
		else if (symbol.equals("Kr"))
			valor = 1.91;
		else if (symbol.equals("Y"))
			valor = 1.78;
		else if (symbol.equals("Zr"))
			valor = 1.56;
		else if (symbol.equals("Nb"))
			valor = 1.48;
		else if (symbol.equals("Ru") || symbol.equals("I"))
			valor = 1.4;
		else if (symbol.equals("Ag"))
			valor = 1.59;
		else if (symbol.equals("Cd"))
			valor = 1.69;
		else if (symbol.equals("In"))
			valor = 1.63;
		else if (symbol.equals("Sn") || symbol.equals("Sb"))
			valor = 1.46;
		else if (symbol.equals("Xe"))
			valor = 1.98;
		else if (symbol.equals("Cs"))
			valor = 1.67;
		else if (symbol.equals("La"))
			valor = 1.87;
		else if (symbol.equals("Ce"))
			valor = 1.83;
		else if (symbol.equals("Pr"))
			valor = 1.82;
		else if (symbol.equals("Nd"))
			valor = 1.81;
		else if (symbol.equals("Pm") || symbol.equals("Sm"))
			valor = 1.8;
		else if (symbol.equals("Eu"))
			valor = 1.99;
		else if (symbol.equals("Gd"))
			valor = 1.79;
		else if (symbol.equals("Tb"))
			valor = 1.76;
		else if (symbol.equals("Dy"))
			valor = 1.75;
		else if (symbol.equals("Ho"))
			valor = 1.74;
		else if (symbol.equals("Er"))
			valor = 1.73;
		else if (symbol.equals("Tm") || symbol.equals("Lu"))
			valor = 1.72;
		else if (symbol.equals("Yb"))
			valor = 1.94;
		else if (symbol.equals("Hf"))
			valor = 1.57;
		else if (symbol.equals("Ta"))
			valor = 1.43;
		else if (symbol.equals("W") || symbol.equals("Os"))
			valor = 1.37;
		else if (symbol.equals("Ir"))
			valor = 1.32;
		else if (symbol.equals("Hg"))
			valor = 1.7;
		else if (symbol.equals("Tl"))
			valor = 1.55;
		else if (symbol.equals("Pb") || symbol.equals("Bi"))
			valor = 1.54;
		else if (symbol.equals("Po"))
			valor = 1.68;
		else if (symbol.equals("At"))
			valor = 1.7;
		else if (symbol.equals("Rn"))
			valor = 2.4;
		else if (symbol.equals("F"))
			valor = 2.0;
		else if (symbol.equals("Ta"))
			valor = 1.9;
		else if (symbol.equals("Ac"))
			valor = 1.88;
		else if (symbol.equals("Th"))
			valor = 1.79;
		else if (symbol.equals("Pa"))
			valor = 1.61;
		else if (symbol.equals("U"))
			valor = 1.58;
		else if (symbol.equals("Np"))
			valor = 1.55;
		else if (symbol.equals("Pu"))
			valor = 1.53;
		else if (symbol.equals("Am"))
			valor = 1.51;
		else if (symbol.equals("Rf") || symbol.equals("Db") || symbol.equals("Sg") || symbol.equals("Bh") || symbol.equals("Hs") || symbol.equals("Mt"))
			valor = 1.6;
		return valor;
	}

}
