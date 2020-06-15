package domain;

import java.util.ArrayList;
import java.util.Arrays;

import javax.vecmath.Point3d;

import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.qsar.IHybridDescriptor.Descriptor;

public class CentroDescriptor {

	IAtomContainer centroDescriptor;
	private TipoCentroDescriptor tipo;
	Secuencia secuencia;
	Point3d centroMasa;
	int id;
	double[] valoresTopologTotal, valoresTopogTotal, valoresTopololTotal, valoresTopolTotal;

	CentroDescriptor(IAtomContainer centroDescriptor, TipoCentroDescriptor tipo) {

		this.centroDescriptor = centroDescriptor;
		this.tipo = tipo;
		id = 0;
		valoresTopogTotal = new double[3];
		valoresTopologTotal = new double[3];
		valoresTopolTotal = new double[3];
		valoresTopololTotal = new double[3];
		Arrays.fill(valoresTopogTotal, Double.MAX_VALUE);
		Arrays.fill(valoresTopologTotal, Double.MAX_VALUE);
		Arrays.fill(valoresTopolTotal, Double.MAX_VALUE);
		Arrays.fill(valoresTopololTotal, Double.MAX_VALUE);
		calcularCentroMasa();
		secuencia = new Secuencia();
	}

	private void calcularCentroMasa() {

		centroMasa = centroDescriptor.getAtom(0).getPoint3d();
		String simb = centroDescriptor.getAtom(0).getSymbol();
		if (!simb.equals("F") && !simb.equals("Cl") && !simb.equals("Br") && !simb.equals("I")) {
			centroMasa = GeometryTools.get3DCentreOfMass(centroDescriptor);
		}
	}

	public void AllIndiceTotalTopographical() {
		if (tipo == TipoCentroDescriptor.ANILLO || tipo == TipoCentroDescriptor.CLUSTER3 || tipo == TipoCentroDescriptor.CLUSTER4)
			AllTotalIndicesTopographicalMultiAtomicos();
		else
			AllTotalIndicesTopographicalMonoAtomicos();
	}
	
	public void AllIndiceTotalTopological() {
		if (tipo == TipoCentroDescriptor.ANILLO || tipo == TipoCentroDescriptor.CLUSTER3 || tipo == TipoCentroDescriptor.CLUSTER4)
			AllTotalIndicesTopologicalMultiAtomicos();
		else
			AllTotalIndicesTopologicalMonoAtomicos();
	}

	private void AllTotalIndicesTopographicalMonoAtomicos() {
		IAtom at = centroDescriptor.getFirstAtom();
		valoresTopogTotal[0] = (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
		valoresTopogTotal[1] = (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
		valoresTopogTotal[2] = (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
	}
	
	private void AllTotalIndicesTopologicalMonoAtomicos() {
		IAtom at = centroDescriptor.getFirstAtom();
		valoresTopolTotal[0] = (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
		valoresTopolTotal[1] = (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
		valoresTopolTotal[2] = (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
	}

	private void AllTotalIndicesTopographicalMultiAtomicos() {
		double electro = 0.0, refracto = 0.0, lipo = 0.0;

		for (IAtom at : centroDescriptor.atoms()) {
			if (!at.getSymbol().equals("H")) {
				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
				refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
				lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
			}
		}
		valoresTopogTotal[0] = electro;
		valoresTopogTotal[1] = refracto;
		valoresTopogTotal[2] = lipo;
	}
	
	private void AllTotalIndicesTopologicalMultiAtomicos() {
		double electro = 0.0, refracto = 0.0, lipo = 0.0;

		for (IAtom at : centroDescriptor.atoms()) {
			if (!at.getSymbol().equals("H")) {
				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
				refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
				lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
			}
		}
		valoresTopolTotal[0] = electro;
		valoresTopolTotal[1] = refracto;
		valoresTopolTotal[2] = lipo;
	}

	private int[] numeroAtomos() {

		int[] result = new int[2];
		Arrays.fill(result, 0);
		for (IAtom a : centroDescriptor.atoms()) {
			if (!a.getSymbol().equals("H")) {
				result[0]++;
				if (!a.getSymbol().equals("C"))
					result[1]++;
			}
		}

		return result;
	}

	public TipoCentroDescriptor getTipo() {
		return tipo;
	}

	public final Point3d getCentroMasa() {
		return centroMasa;
	}

	public final int getId() {
		return id;
	}

	public final void setId(int id) {
		this.id = id;
	}

	public final String getSecuencia() {
		return secuencia.getCodigo(centroDescriptor, tipo);
	}

	public final String getNombreCD() {
		int[] atomos = numeroAtomos();
		return (tipo.toString().equals("ANILLO")) ? tipo.toString() + atomos[0] : tipo.toString();
	}

	public double[] getValoresTopogAtomo(int pos) {
		double[] valor = { 0, 0, 0 };

		valor[0] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.ELECTROTOPOGRAPHIC);
		valor[1] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.REFRACTEDTOPOGRAPHIC);
		valor[2] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.LIPOTOPOGRAPHIC);

		return valor;
	}
	
	public double[] getValoresTopolAtomo(int pos) {
		double[] valor = { 0, 0, 0 };

		valor[0] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.ELECTROTOPOLOGICAL);
		valor[1] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.REFRACTEDTOPOLOGICAL);
		valor[2] = (Double) centroDescriptor.getAtom(pos).getProperty(Descriptor.LIPOTOPOLOGICAL);

		return valor;
	}

	public int getCantidadAtomos() {
		return centroDescriptor.getAtomCount();
	}
	
	public Iterable<IAtom> getListAtom() {
		return centroDescriptor.atoms();
	}

	public String getSimbolAtomo(int pos) {
		return centroDescriptor.getAtom(pos).getSymbol();
	}

	public String getIdAtomo(int pos) {
		return centroDescriptor.getAtom(pos).getID();
	}

	public double[] getValoresTopogTotal(char indice) {
		switch (indice) {
			case 'E':
				return new double[] { valoresTopogTotal[0] };
			case 'R':
				return new double[] { valoresTopogTotal[1] };
			case 'L':
				return new double[] { valoresTopogTotal[2] };
			default:
				return valoresTopogTotal;
		}

	}
	
	public double getValorTotalVector(String index){
		double valor =0;
		if(index.equals("E"))
			valor = valoresTopogTotal[0];
		else if(index.equals("R")) 
			valor = valoresTopogTotal[1];
		else if(index.equals("L"))
			valor = valoresTopogTotal[2];
		else if(index.equals("ER"))
			valor = Math.sqrt(Math.pow(valoresTopogTotal[0], 2) + Math.pow(valoresTopogTotal[1], 2));
		else if(index.equals("EL"))
			valor = Math.sqrt(Math.pow(valoresTopogTotal[0], 2) + Math.pow(valoresTopogTotal[2], 2));
		else if(index.equals("RL"))
			valor = Math.sqrt(Math.pow(valoresTopogTotal[1], 2) + Math.pow(valoresTopogTotal[2], 2));
		else if(index.equals("ERL"))
			valor = Math.sqrt(Math.pow(valoresTopogTotal[0], 2) + Math.pow(valoresTopogTotal[1], 2) + Math.pow(valoresTopogTotal[2], 2));
		return valor;
	}
	
	public double[] getValoresTopogTotalNorm(char indice) {
		double electro = 0.0, refracto = 0.0, lipo = 0.0;
		double[] valoresTopolTotalNorm = new double[3];
		for (IAtom at : centroDescriptor.atoms()) {
			if (!at.getSymbol().equals("H")) {
				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM);
				refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM);
				lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM);
			}
		}
		valoresTopolTotalNorm[0] = electro;
		valoresTopolTotalNorm[1] = refracto;
		valoresTopolTotalNorm[2] = lipo;
		
		switch (indice) {
			case 'E':
				return new double[] { valoresTopolTotalNorm[0] };
			case 'R':
				return new double[] { valoresTopolTotalNorm[1] };
			case 'L':{
				return new double[] { valoresTopolTotalNorm[2] };
			}
			default:{
				return new double[] {valoresTopolTotalNorm[0], valoresTopolTotalNorm[1], valoresTopolTotalNorm[2]};
			}
		}

	}
	
	public double[] getValoresTopolTotal(char indice) {
		switch (indice) {
			case 'E':
				return new double[] { valoresTopolTotal[0] };
			case 'R':
				return new double[] { valoresTopolTotal[1] };
			case 'L':
				return new double[] { valoresTopolTotal[2] };
			default:
				return valoresTopolTotal;
		}

	}

	public double[] getValoresTopogTotal() {
		return valoresTopogTotal;
	}
	
	public double[] getValoresTopolTotal() {
		return valoresTopolTotal;
	}

	// *****************************************************************************************************
	// ************************IMPLEMENTADO POR JUAN LUIS PANEQUE*******************************************
	// ******************************************************************************************************

	public String getIdentifier() {
		if (tipo == TipoCentroDescriptor.ANILLO)
			return getTipo().toString() + getCantidadAtomos();
		else if (tipo == TipoCentroDescriptor.CLUSTER3 || tipo == TipoCentroDescriptor.CLUSTER4)
			return getTipo().toString();
		else
			return getTipo().toString() + "-" + getSecuencia();
	}
	
	public ArrayList<IAtom> listAtomsDifferent(CentroDescriptor cd){
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		for(IAtom at : centroDescriptor.atoms())
			if(!at.getSymbol().equals("H"))
				listAtom.add(at);
		
		for(IAtom at : cd.getListAtom()){
			if(!at.getSymbol().equals("H")){
				boolean flag = false;
				String atomo1 = at.getSymbol()+at.getID();
				int j = 0;
				while(!flag && j<listAtom.size()){
					if(!listAtom.get(j).getSymbol().equals("H")){
						String atomo2 = listAtom.get(j).getSymbol()+listAtom.get(j).getID();
						if(atomo1.equals(atomo2))
							flag = true; 
					}
					j++;	
				}
				if(!flag)
					listAtom.add(at);
			}	
		}
		return listAtom;
	}

}
