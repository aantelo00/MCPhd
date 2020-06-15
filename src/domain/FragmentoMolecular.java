package domain;

import java.util.ArrayList;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.qsar.IHybridDescriptor.Descriptor;

public class FragmentoMolecular {

	private ArrayList<CentroDescriptor> centDescriptores;
	private ArrayList<IAtom> centAtoms;
    private String tipo, tipoAtoms;
	private int gradoFragmento;
	private double indexSimilary;
	private Molecule molecula;

	public FragmentoMolecular() {
		centDescriptores = new ArrayList<CentroDescriptor>();
		centAtoms = new ArrayList<IAtom>();
		gradoFragmento = 0;
		tipo = "";
		tipoAtoms = "";
		//indexSimilary = 0.0;
		molecula = new Molecule();
	}

	/**
	 * @return the molecula
	 */
	public Molecule getMolecula() {
		return molecula;
	}

	/**
	 * @param molecula
	 *            the molecula to set
	 */
	public void setMolecula(Molecule molecula) {
		this.molecula = molecula;
	}
	
	public void setTipo(String tipo) {
		this.tipo = tipo;
	}
	
	public void setGradoFragmento(int gradoFragmento) {
		this.gradoFragmento = gradoFragmento;
	}

	@Override
	public boolean equals(Object arg0) {
		FragmentoMolecular fragment = (FragmentoMolecular) arg0;
		return (molecula.getNombre().equals(fragment.getMolecula().getNombre())) && (isSecuenceCorrect(fragment.getSecuenciaFragmentoCD()));
	}

	@Override
	public int hashCode() {
		//return getSecuenciaFragmento().hashCode();
		return getTipo().hashCode();
	}

	public double[] getDistancia(char method) {
		if (gradoFragmento == 0)
			return new double[] { 0.0 };

		double[] distancia = null;
		switch (method) {
			case 'L': {
				int n = gradoFragmento - 1;
				distancia = new double[n];
				int k = 0;
				for (int i = 0 ; i < n; i++) {
					if(!(i == n-1))
					  distancia[k] = centDescriptores.get(i).getCentroMasa().distance(centDescriptores.get(i+1).getCentroMasa());
					else
						distancia[k] = centDescriptores.get(i).getCentroMasa().distance(centDescriptores.get(0).getCentroMasa());	
					k++;
				}
				break;
			}
			case 'T': {
				int n = gradoFragmento * (gradoFragmento - 1) / 2;
				distancia = new double[n];

				for (int i = 0 ; i < gradoFragmento - 1; i++) {
					for (int j = i + 1; j < gradoFragmento; j++) {
						distancia[i] = centDescriptores.get(i).getCentroMasa().distance(centDescriptores.get(j).getCentroMasa());
					}
				}

				break;
			}
			case 'C': {
				int n = gradoFragmento;
				distancia = new double[n];
				for (int i = 0 ; i < gradoFragmento; i++) {
					if (i == gradoFragmento - 1) {
						distancia[i] = centDescriptores.get(i).getCentroMasa().distance(centDescriptores.get(0).getCentroMasa());

					} else {
						distancia[i] = centDescriptores.get(i).getCentroMasa().distance(centDescriptores.get(i + 1).getCentroMasa());
					}
				}

				break;
			}
		}
		return distancia;
	}
	
	public void addAtom(IAtom atom){
		if(!atom.getSymbol().equals("H")){
			boolean flag = true;
			String atom1 = atom.getSymbol() + atom.getID();
			for(IAtom at : centAtoms){
				String atom2 = at.getSymbol() + at.getID();
				if(atom1.equals(atom2))
					flag = false;
			}
			if(flag){
				if (centAtoms.size() > 0)
					   tipoAtoms += ",";
				    centAtoms.add(atom);
				    tipoAtoms += atom.getSymbol() + atom.getID();
			}   
		} 
	}

	public void addCentroDescriptor(CentroDescriptor centro) {
		centDescriptores.add(centro);
		for(IAtom at : centro.getListAtom()){
			addAtom(at);
		}
		
		gradoFragmento++;

		if (gradoFragmento > 1)
			tipo += "-";

		if (centro.getTipo().toString().equals("ANILLO") || centro.getTipo().toString().equals("CLUSTER"))
			tipo += centro.getTipo().toString() + centro.getCantidadAtomos() + "[" + centro.getSecuencia() + "]";
		else
			tipo += centro.getTipo().toString() + "[" + centro.getSecuencia() + "]";
	}
	
	public void addAtoms(IAtom atom) {
		if (centAtoms.size() > 0)
			tipoAtoms += ",";
		centAtoms.add(atom);
		tipoAtoms += atom.getSymbol() + atom.getID();
	}	

	public String getSecuenciaFragmentoCD() {
		String secuencia = centDescriptores.get(0).getSecuencia();
		for (int i = 1; i < gradoFragmento; i++) {
			secuencia += "," + centDescriptores.get(i).getSecuencia();
		}
		return secuencia;
	}
	
	public String getSecuenciaFragmentoAtom() {
		String secuencia = centAtoms.get(0).getSymbol() + centAtoms.get(0).getID() ;
		for (int i = 1; i < getCantAtoms(); i++) {
			secuencia += "," + centAtoms.get(i).getSymbol() + centAtoms.get(i).getID() ;
		}
		return secuencia;
	}

	public boolean isSecuenceCorrect(String sec) {
		String[] array = sec.split(",");
		String secuencia = getSecuenciaFragmentoCD();
		for (String s : array) {
			if (secuencia.indexOf(s) == -1)
				return false;
		}
		return true;
	}

	public void setIndexSimilary(double indexSimilary) {
		this.indexSimilary = indexSimilary;
	}

	public ArrayList<CentroDescriptor> getCentDescriptores() {
		return centDescriptores;
	}

	public String getTipo() {
		return tipo;
	}
	
	public String getTipoAtoms() {
		return tipoAtoms;
	}

	public int getGradoFragmento() {
		return gradoFragmento;
	}
	
	public int getCantAtoms() {
		return centAtoms.size();
	}

	public double getIndexSimilary() {
		return indexSimilary;
	}

	public double[] getValoresTopograficosTotal() {
		double[] valores = {0.0,0.0,0.0,0.0};
		for (int i = 0; i < gradoFragmento; i++) {
			//System.out.println(centDescriptores.get(i).getSecuencia());
			valores[0] += centDescriptores.get(i).getValoresTopogTotal()[0];
			valores[1] += centDescriptores.get(i).getValoresTopogTotal()[1];
			valores[2] += centDescriptores.get(i).getValoresTopogTotal()[2];
		}
		valores[3] = getDistanciaTotal();
		return valores;
	}
	
	private ArrayList<IAtom> listAtomsCD(){
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		if(gradoFragmento > 1){
			for (int i = 0; i < gradoFragmento-1; i++)
				   for(int j=i+1; j< gradoFragmento; j++)
					  listAtom.addAll(centDescriptores.get(i).listAtomsDifferent(centDescriptores.get(j)));
		} else{
			for(IAtom a : centDescriptores.get(0).getListAtom())
				listAtom.add(a);
		}
		
		return listAtom;
	}
	
	private ArrayList<IAtom> listAtomsAtom(){
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		for (int i = 0; i < getCantAtoms(); i++)
			listAtom.add(centAtoms.get(i));
		return listAtom;
	}
	
	public double getValoresTopograficosPMCCD(String index){
		double valor =0.0;
		ArrayList<IAtom> listAtom = listAtomsCD();
		ArrayList<IAtom> listAtoms = new ArrayList<IAtom>();
		
		for(IAtom at : listAtom){
			boolean flag = false;
			String atomo1 = at.getSymbol() + at.getID();
			int j = 0;
			while(!flag && j<listAtoms.size()){
			    String atomo2 = listAtoms.get(j).getSymbol() + listAtoms.get(j).getID();
			    if(atomo1.equals(atomo2))
					flag = true; 
				j++;	
			}
			if(!flag)
				listAtoms.add(at);
	    }
		
		for(int i=0; i<listAtoms.size(); i++)
			if(!listAtom.get(i).getSymbol().equals("H"))
				if(index == "E")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM);
				else if(index == "R")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM);
				else if(index == "L")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM);
				else if(index == "ER") 
					   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
							    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2));
				else if(index == "EL") 
				   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
				else if(index == "RL")
					   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) +
							    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
				else
				   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) + 
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
		
		    return valor;
	}
	
	public double getValoresTopograficosPMCAtoms(String index){
		double valor =0.0;
		ArrayList<IAtom> listAtom = listAtomsAtom();
		ArrayList<IAtom> listAtoms = new ArrayList<IAtom>();
		
		for(IAtom at : listAtom){
			boolean flag = false;
			String atomo1 = at.getSymbol() + at.getID();
			int j = 0;
			while(!flag && j<listAtoms.size()){
			    String atomo2 = listAtoms.get(j).getSymbol() + listAtoms.get(j).getID();
			    if(atomo1.equals(atomo2))
					flag = true; 
				j++;	
			}
			if(!flag)
				listAtoms.add(at);
	    }
		
		for(int i=0; i<listAtoms.size(); i++)
			if(!listAtom.get(i).getSymbol().equals("H"))
				if(index == "E")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM);
				else if(index == "R")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM);
				else if(index == "L")
				   valor += (Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM);
				else if(index == "ER") 
					   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
							    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2));
				else if(index == "EL") 
				   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
				else if(index == "RL")
					   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) +
							    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
				else
				   valor += Math.sqrt(Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) + 
						    Math.pow((Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHICNORM), 2));
		
		    return valor;
	}
	
	
	
	public double[] getValoresTopograficosPMC() {
		double[] valores = {0.0,0.0,0.0};
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		ArrayList<IAtom> listAtoms = new ArrayList<IAtom>();
		if(gradoFragmento > 1){
			for (int i = 0; i < gradoFragmento-1; i++)
				   for(int j=i+1; j< gradoFragmento; j++)
					  listAtom.addAll(centDescriptores.get(i).listAtomsDifferent(centDescriptores.get(j)));
		} else{
			for(IAtom a : centDescriptores.get(0).getListAtom())
				listAtom.add(a);
		}
			
		
	
		for(IAtom at : listAtom){
				boolean flag = false;
				String atomo1 = at.getSymbol() + at.getID();
				int j = 0;
				while(!flag && j<listAtoms.size()){
				    String atomo2 = listAtoms.get(j).getSymbol() + listAtoms.get(j).getID();
				    if(atomo1.equals(atomo2))
						flag = true; 
					j++;	
				}
				if(!flag)
					listAtoms.add(at);
		}
		
		@SuppressWarnings("unused")
		String cadena = "";
		for(IAtom at : listAtoms)
			cadena += at.getSymbol() + at.getID() + ",";
		
		//System.out.println(cadena);
		
		for(int i=0; i<listAtoms.size(); i++){
			if(!listAtom.get(i).getSymbol().equals("H")){
				valores[0] += (Double) listAtoms.get(i).getProperty(Descriptor.ELECTROTOPOGRAPHIC);
				valores[1] += (Double) listAtoms.get(i).getProperty(Descriptor.REFRACTEDTOPOGRAPHICNORM);
				valores[2] += (Double) listAtoms.get(i).getProperty(Descriptor.LIPOTOPOGRAPHIC);
			}
			
		}
		
		return valores;
	}	

	public double[] vectorValoresTopograficos() {
		double[] vector = new double[gradoFragmento * 3 + 1];
		int c = 0;
		for (int i = 0; i < gradoFragmento; i++) {
			vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[0];
			vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[1];
			vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[2];
		}
		vector[c] = getDistanciaTotal();
		return vector;
	}
	
	public double[] vectorValoresTopologicos() {
		double[] vector = new double[gradoFragmento * 3 + 1];
		int c = 0;
		for (int i = 0; i < gradoFragmento; i++) {
			vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[0];
			vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[1];
			vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[2];
		}
		vector[c] = getDistanciaTotal();
		return vector;
	}

	private double getDistanciaTotal() {
		if (gradoFragmento == 1) {
			return 0;
		} else {
			double[] distancias = getDistancia('T');
			int n = gradoFragmento * (gradoFragmento - 1) / 2;
			double dt = 0.0;
			for (double d : distancias)
				dt += d;
			return dt / n;
		}
	}

	public double[] vectorValoresTopograficos(char indice) {
		double[] vector = null;
		int c = 0;
		switch (indice) {
			case 'A': {
				vector = new double[gradoFragmento * 3 + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[0];
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[1];
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[2];
				}
				break;
			}
			case 'E': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[0];
				}
				break;
			}
			case 'R': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[1];
				}
				break;
			}
			case 'L': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopogTotal()[2];
				}
				break;
			}
		}
		vector[0] = centDescriptores.get(0).getCentroMasa().distanceSquared(centDescriptores.get(1).centroMasa);
		return vector;
	}
	
	public double[] vectorValoresTopologicos(char indice) {
		double[] vector = null;
		int c = 0;
		switch (indice) {
			case 'A': {
				vector = new double[gradoFragmento * 3 + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[0];
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[1];
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[2];
				}
				break;
			}
			case 'E': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[0];
				}
				break;
			}
			case 'R': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[1];
				}
				break;
			}
			case 'L': {
				vector = new double[gradoFragmento + 1];
				for (int i = 0; i < gradoFragmento; i++) {
					vector[c++] = centDescriptores.get(i).getValoresTopolTotal()[2];
				}
				break;
			}
		}
		vector[c] = getDistanciaTotal();
		return vector;
	}

	public double getSimilarityOrder2(TypeSimilaryFunction funcSim, boolean max, FragmentoMolecular fragment, char indice) throws CDKException {
		double[] frag1 = vectorValoresTopograficos(indice);
		double[] frag2 = fragment.vectorValoresTopograficos(indice);

		double index1 = new SimilaryFunctions().CalculeSimilaryFunction(funcSim, frag1, frag2);

		double[] frag3;
		if (frag2.length == 7)
			frag3 = new double[] { frag2[3], frag2[4], frag2[5], frag2[0], frag2[1], frag2[2], frag2[6] };
		else
			frag3 = new double[] { frag2[1], frag2[0], frag2[2] };

		double index2 = new SimilaryFunctions().CalculeSimilaryFunction(funcSim, frag1, frag3);

		if (max)
			return Math.max(index1, index2);
		else
			return Math.min(index1, index2);
	}
	
	public ArrayList<IAtom> getListAtomsCD(){
		return centAtoms;
	}
	
	public ArrayList<IAtom> getListAtomsAtom(){
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		for(int i=0; i< centAtoms.size(); i++){
				listAtom.add(centAtoms.get(i));
		}
		return listAtom;
	}
	
	public boolean ContainCenterDescriptor(CentroDescriptor center){
		boolean flag = false;
		for(int i=0; i<centDescriptores.size(); i++){
			if(centDescriptores.get(i).getSecuencia().equals(center.getSecuencia()))
				flag=true;
		}
		return flag;
	}
	
	public boolean ContainAtom(IAtom center){
		boolean flag = false;
		for(int i=0; i<centAtoms.size(); i++){
			String atom1 = center.getSymbol() + center.getID();
			String atom2 = centAtoms.get(i).getSymbol() + centAtoms.get(i).getID();
			if(atom1.equals(atom2))
				flag=true;
		}
		return flag;
	}
}
