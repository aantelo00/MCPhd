package domain;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import javax.vecmath.Point3d;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.ringsearch.SSSRFinder;

public class ContenedorCD {

	private IAtomContainer molecula;
	private ArrayList<IAtom> metilo = new ArrayList<IAtom>();
	private ArrayList<IAtom> metileno = new ArrayList<IAtom>();
	private ArrayList<IAtom> metino = new ArrayList<IAtom>();
	private ArrayList<IAtom> clusterOrden3 = new ArrayList<IAtom>();
	private ArrayList<IAtom> clusterOrden4 = new ArrayList<IAtom>();
	private ArrayList<IAtom> heteroAtomo = new ArrayList<IAtom>();
	private IAtomContainerSet metiloSet, metilenoSet, metinoSet, clusterOrden3Set,
	clusterOrden4Set, heteroAtomoSet;
	private TreeMap<CDescrip, IAtomContainerSet> cdescriptores;
	private int ID;
	private double[][] distanciasTopograficas;
	//ArrayList<List<List<IAtom>>> caminos;
	private TreeMap <String, List<List<IAtom>>> caminos;


	public enum CDescrip{
		ANILLO,
		CLUSTER4,
		CLUSTER3,
		METINO,
		METILENO,
		METILO,
		HETEROATOMO		
	}


	public ContenedorCD(IAtomContainer m) {

		this.molecula = m;
		cdescriptores = new TreeMap<CDescrip,IAtomContainerSet>();
		//caminos = new ArrayList<List<List<IAtom>>>();
		caminos = new TreeMap<String, List<List<IAtom>>>();
		ID = 0;


		// Faltan los caminos

	}

	public ContenedorCD() {

		cdescriptores = new TreeMap<CDescrip,IAtomContainerSet>();
		//caminos = new ArrayList<List<List<IAtom>>>();
		caminos = new TreeMap<String, List<List<IAtom>>>();
		ID = 0;
		// Faltan los caminos

	}

//	public void testAtomContainer() {
//		IAtomContainer partContainer = DefaultChemObjectBuilder.getInstance().newAtomContainer();
//		partContainer.addAtom(clusterOrden3.get(0));
//		List<IBond> parts = this.molecula.getConnectedBondsList(clusterOrden3.get(0));
//		for (IBond aBond : parts) {
//			Iterator<IAtom> bondedAtoms = aBond.atoms().iterator();
//			while (bondedAtoms.hasNext()) {
//				IAtom bondedAtom = bondedAtoms.next();
//				if (!partContainer.contains(bondedAtom)
//						&& !bondedAtom.getSymbol().equals("H"))
//					partContainer.addAtom(bondedAtom);
//			}
//			partContainer.addBond(aBond);
//		}
//
//		IAtom atomoIn = partContainer.getAtom(0);
//		System.out.print(atomoIn.getSymbol() + atomoIn.getID());
//		for (int i = 1; i < partContainer.getAtomCount(); i++) {
//			IAtom at = partContainer.getAtom(i);
//			IBond bo = partContainer.getBond(atomoIn, at);
//			String orden = "";
//			switch (bo.getOrder()) {
//			case SINGLE:
//				orden = "s";
//				break;
//			case DOUBLE:
//				orden = "d";
//				break;
//			case TRIPLE:
//				orden = "t";
//				break;
//
//			}
//			System.out.print(orden + at.getSymbol() + at.getID());
//		}
//		System.out.println("");
//
//		// try {
//		// List<String> types = new ArrayList<String>();
//		// types.add("hybrid");
//		// DescriptorEngine engine = new DescriptorEngine(types);
//		// engine.process(molecula);
//		// } catch (CDKException e) {
//		// // TODO Auto-generated catch block
//		// e.printStackTrace();
//		// }
//
//	}

	private void updateTotalIndicesTopologicos(IAtomContainerSet containers) {

		double electro = 0.0, refracto = 0.0, lipo = 0.0;
		for (IAtomContainer cont : containers.atomContainers()) {
			electro = refracto = lipo = 0.0;
			for (IAtom at : cont.atoms()) {
				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
				refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
				lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
			}

			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL,electro);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL,refracto);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL,lipo);
		}

	}

	private void updateTotalIndicesTopograficos(IAtomContainerSet containers) {

		double electro = 0.0, refracto = 0.0, lipo = 0.0;
		for (IAtomContainer cont : containers.atomContainers()) {
			electro = refracto = lipo = 0.0;
			for (IAtom at : cont.atoms()) {
				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
				refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
				lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
			}

			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC,electro);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC,refracto);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC,lipo);
		}

	}

	private void updateTotalIndicesTopologicosMonoAtomicos(IAtomContainerSet containers) {

		double electro = 0.0, refracto = 0.0, lipo = 0.0;
		for (IAtomContainer cont : containers.atomContainers()) {
			electro = refracto = lipo = 0.0;
			IAtom at = cont.getFirstAtom();
			//for (IAtom at : cont.atoms()) {
			electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
			refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
			lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
			//}

			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL,electro);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL,refracto);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL,lipo);
		}

	}

	private void updateTotalIndicesTopograficosMonoAtomicos(IAtomContainerSet containers) {

		double electro = 0.0, refracto = 0.0, lipo = 0.0;
		for (IAtomContainer cont : containers.atomContainers()) {
			electro = refracto = lipo = 0.0;
			IAtom at = cont.getFirstAtom();
			//for (IAtom at : cont.atoms()) {
			electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
			refracto += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
			lipo += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
			//}

			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC,electro);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC,refracto);
			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC,lipo);
		}

	}

//	private void updateTotalIndices(IAtomContainerSet containers) {
//
//		double electro = 0.0, refracto = 0.0, lipo = 0.0;
//		double electro1 = 0.0, refracto1 = 0.0, lipo1 = 0.0;
//		for (IAtomContainer cont : containers.atomContainers()) {
//			electro = refracto = lipo = 0.0;
//			electro1 = refracto1 = lipo1 = 0.0;
//			for (IAtom at : cont.atoms()) {
//				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
//				refracto += (Double) at	.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
//				electro += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
//				electro1 += (Double) at.getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
//				refracto1 += (Double) at.getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
//				lipo1 += (Double) at.getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
//
//			}
//
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL, electro);
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL,refracto);
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL,lipo);
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC,electro1);
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC,refracto1);
//			cont.setProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC,lipo1);
//		}
//
//	}



	public void indiceTotalTopologico(ContenedorCD.CDescrip tipo){

		//    	switch(tipo){
		//    	case HETEROATOMO: updateTotalIndicesTopologicos(heteroAtomoSet); break;
		//    	case CLUSTER3: updateTotalIndicesTopologicos(clusterOrden3Set); break;
		//    	case CLUSTER4: updateTotalIndicesTopologicos(clusterOrden4Set); break;
		//    	case METILO:  updateTotalIndicesTopologicos(metiloSet); break;
		//    	case METINO: updateTotalIndicesTopologicos(metinoSet); break;
		//    	case METILENO: updateTotalIndicesTopologicos(metilenoSet); break;
		//    	case ANILLO: updateTotalIndicesTopologicos(anillos); break;
		//
		//    	}
		if(tipo == CDescrip.ANILLO || tipo == CDescrip.CLUSTER3 || tipo == CDescrip.CLUSTER4)
			updateTotalIndicesTopologicos(cdescriptores.get(tipo));
		else 
			updateTotalIndicesTopologicosMonoAtomicos(cdescriptores.get(tipo));		
	} 

	public void indiceTotalTopografico(ContenedorCD.CDescrip tipo){

		//		switch(tipo){
		//		case HETEROATOMO: updateTotalIndicesTopograficos(heteroAtomoSet); break;
		//		case CLUSTER3: updateTotalIndicesTopograficos(clusterOrden3Set); break;
		//		case CLUSTER4: updateTotalIndicesTopograficos(clusterOrden4Set); break;
		//		case METILO:  updateTotalIndicesTopograficos(metiloSet); break;
		//		case METINO: updateTotalIndicesTopograficos(metinoSet); break;
		//		case METILENO: updateTotalIndicesTopograficos(metilenoSet); break;
		//		case ANILLO: updateTotalIndicesTopograficos(anillos); break;
		//
		//		}

		if(tipo == CDescrip.ANILLO || tipo == CDescrip.CLUSTER3 || tipo == CDescrip.CLUSTER4)
			updateTotalIndicesTopograficos(cdescriptores.get(tipo)); 
		else 
			updateTotalIndicesTopograficosMonoAtomicos(cdescriptores.get(tipo));		
	} 



	public String codificacionCentroDescriptor(IAtomContainerSet containers, String nombre){
		StringBuffer codigos = new StringBuffer();

		for(IAtomContainer partContainer : containers.atomContainers()){
			String tipo = partContainer.getProperty("Tipo").toString() + "_";
			int n = numeroHeteroAtomos(partContainer);
			int cant = partContainer.getAtomCount();
			String cadena = cant + "_" + n + "_" + tipo;
			codigos.append(nombre + "_" + cadena);
			IAtom atomoIn = partContainer.getAtom(0);
			codigos.append(atomoIn.getSymbol() + atomoIn.getID());
			for (int i = 1; i < partContainer.getAtomCount(); i++) {
				IAtom at = partContainer.getAtom(i);
				IBond bo = partContainer.getBond(atomoIn, at);
				String orden = "";
				if(bo.getFlag(CDKConstants.ISAROMATIC))
					orden = "a";
				else{
					switch (bo.getOrder()) {
					case SINGLE:
						orden = "s";
						break;
					case DOUBLE:
						orden = "d";
						break;
					case TRIPLE:
						orden = "t";
						break;
					default: 
						orden = "unset";
						break;
					}
				}

				codigos.append(orden + at.getSymbol() + at.getID());
			}
			codigos.append("_");
			NumberFormat format = NumberFormat.getInstance();
			format.setMaximumFractionDigits(5);
			double refratedTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL);
			double electroTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL);
			double lipoTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL);
			double refratedTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC);
			double electroTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC);
			double lipoTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC);
			String valores = format.format(refratedTopol) + "_" +  format.format(electroTopol) + "_" + format.format(lipoTopol) + "_"
			+ format.format(refratedTopog) + "_" + format.format(electroTopog) + "_" + format.format(lipoTopog);


			codigos.append(valores + "\n");
			//codigos.append("\n");
		}
		codigos.trimToSize();
		return codigos.toString();

	} 

	public String codificarAnillo(IAtomContainerSet atomContainer, String nombre) {

		StringBuffer stb = new StringBuffer();
		//stb.append(nombre + "_" );
		for(IAtomContainer ac : atomContainer.atomContainers()){
			String tipo = ac.getProperty("Tipo").toString() + "_";
			int n = numeroHeteroAtomos(ac);
			int cant = ac.getAtomCount();
			String cadena = cant + "_" + n + "_" + tipo;
			stb.append(nombre + "_" + cadena);
			IRingSet anillos = new SSSRFinder(ac).findEssentialRings();
			for(IAtomContainer an : anillos.atomContainers()){
				ArrayList<IAtom> vertices = new ArrayList<IAtom>();
				ArrayList<IBond> enlaces = new ArrayList<IBond>();
				IAtom firstAtom = an.getFirstAtom();
				stb.append(firstAtom.getSymbol() + firstAtom.getID());
				dephtFirstSearch(an, firstAtom, vertices,enlaces,stb);

			}

			stb.append(valoresIndicesAtomContainer(ac));
			stb.append("\n");
		}
		stb.trimToSize();
		return stb.toString();
	}

	private void dephtFirstSearch(IAtomContainer atomContainer,IAtom firstAtom, ArrayList<IAtom> vertices,
			ArrayList<IBond> enlaces, StringBuffer stb) {
		vertices.add(firstAtom);

		for(IBond bond : atomContainer.getConnectedBondsList(firstAtom)){
			if(!enlaces.contains(bond)){
				enlaces.add(bond);
				IAtom second = bond.getConnectedAtom(firstAtom);
				String cad = codificarAtomo(bond, second);
				stb.append(cad);

				if(!vertices.contains(second))
					dephtFirstSearch(atomContainer, second, vertices, enlaces, stb);
			}
		}

	}

	//CentroDescriptor
	private String valoresIndicesAtomContainer(IAtomContainer partContainer){
		StringBuffer codigos = new StringBuffer();
		codigos.append("_");    
		NumberFormat format = NumberFormat.getInstance();
		format.setMaximumFractionDigits(5);
		double refratedTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL);
		double electroTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL);
		double lipoTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL);
		double refratedTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC);
		double electroTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC);
		double lipoTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC);
		String valores = format.format(refratedTopol) + "_" +  format.format(electroTopol) + "_" + format.format(lipoTopol) + "_"
		+ format.format(refratedTopog) + "_" + format.format(electroTopog) + "_" + format.format(lipoTopog);
		codigos.append(valores);
		codigos.trimToSize();
		return codigos.toString();
	}

	//	public String codificacionAnillo(IAtomContainerSet containers, String nombre){
	//		StringBuffer codigos = new StringBuffer();
	//
	//
	//		for(IAtomContainer partContainer : containers.atomContainers()){
	//
	//			String tipo = partContainer.getProperty("Tipo").toString() + "_";
	//			int n = numeroHeteroAtomos(partContainer);
	//			int cant = partContainer.getAtomCount();
	//			String cadena = cant + "_" + n + "_" + tipo;
	//			codigos.append(nombre + "_" + cadena);
	//
	//			IAtom atomoInc = partContainer.getAtom(0);
	//
	//			for(IAtom at : partContainer.atoms()){
	//				if(at != atomoInc){
	//					String c = codificarAtomo(partContainer, atomoInc, at);
	//					codigos.append(c);
	//				}
	//				atomoInc = at;
	//
	//
	//			}
	//
	//			//			codigos.append(atomoInc.getSymbol() + atomoInc.getID());
	//			//			for (int i = 1; i < cant; i++) {
	//			//				IAtom atomoIn, at;
	//			//
	//			//				if(i == partContainer.getAtomCount() - 1){
	//			//					atomoIn = partContainer.getAtom(i);
	//			//					at = partContainer.getAtom(0);
	//			//				}
	//			//				else{
	//			//					atomoIn = partContainer.getAtom(i);
	//			//					at = partContainer.getAtom(i + 1);
	//			//				}
	//			//
	//			//				IBond bo = partContainer.getBond(atomoIn, at);
	//			//				String orden = "";
	//			//				if(bo.getFlag(CDKConstants.ISAROMATIC))
	//			//					orden = "a";
	//			//				else{
	//			//					switch (bo.getOrder()) {
	//			//					case SINGLE:
	//			//						orden = "s";
	//			//						break;
	//			//					case DOUBLE:
	//			//						orden = "d";
	//			//						break;
	//			//					case TRIPLE:
	//			//						orden = "t";
	//			//						break;
	//			//					default: 
	//			//						orden = "unset";
	//			//						break;
	//			//					}
	//			//				}
	//
	//			//codigos.append(orden + at.getSymbol() + at.getID());
	//			//}
	//			codigos.append("_");    
	//			NumberFormat format = NumberFormat.getInstance();
	//			format.setMaximumFractionDigits(5);
	//			double refratedTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOLOGICAL);
	//			double electroTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOLOGICAL);
	//			double lipoTopol = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOLOGICAL);
	//			double refratedTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALREFRACTEDTOPOGRAPHIC);
	//			double electroTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC);
	//			double lipoTopog = (Double)partContainer.getProperty(IHybridDescriptor.Descriptor.TOTALLIPOTOPOGRAPHIC);
	//			String valores = format.format(refratedTopol) + "_" +  format.format(electroTopol) + "_" + format.format(lipoTopol) + "_"
	//			+ format.format(refratedTopog) + "_" + format.format(electroTopog) + "_" + format.format(lipoTopog);
	//			codigos.append(valores + "\n");
	//		}
	//		codigos.trimToSize();
	//		return codigos.toString();
	//
	//	}

	private String codificarAtomo(IBond bond, IAtom at2){

		//String atomos = at1.getSymbol() + at1.getID() + at2.getSymbol() + at2.getID(); 
		String orden = "";
		if(bond.getFlag(CDKConstants.ISAROMATIC))
			orden = "a";
		else{
			switch (bond.getOrder()) {
			case SINGLE:
				orden = "s";
				break;
			case DOUBLE:
				orden = "d";
				break;
			case TRIPLE:
				orden = "t";
				break;
			default: 
				orden = "unset";
				break;
			}
		}

		String codigos = orden + at2.getSymbol() + at2.getID();
		return codigos;
	}



	//	private String codificarAtomo(IAtomContainer ac,  IAtom at1, IAtom at2){
	//		IBond bo = ac.getBond(at1, at2);
	//		String atomos = at1.getSymbol() + at1.getID() + at2.getSymbol() + at2.getID(); 
	//		String orden = "";
	//		if(bo.getFlag(CDKConstants.ISAROMATIC))
	//			orden = "a";
	//		else{
	//			switch (bo.getOrder()) {
	//			case SINGLE:
	//				orden = "s";
	//				break;
	//			case DOUBLE:
	//				orden = "d";
	//				break;
	//			case TRIPLE:
	//				orden = "t";
	//				break;
	//			default: 
	//				orden = "unset";
	//				break;
	//			}
	//		}
	//
	//		String codigos = at1.getSymbol() + at2.getID() + orden + at2.getSymbol() + at2.getID();
	//		return codigos;
	//	}


	private IAtomContainerSet setVecinos(ArrayList<IAtom> atomocentral, CDescrip tipo) {

		IAtomContainerSet centrosDescriptores = new AtomContainerSet();

		for (int i = 0; i < atomocentral.size(); i++) {
			IAtomContainer partContainer = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
			partContainer.addAtom(atomocentral.get(i));
			List<IBond> parts = this.molecula.getConnectedBondsList(atomocentral.get(i));
			for (IBond aBond : parts) {
				Iterator<IAtom> bondedAtoms = aBond.atoms().iterator();
				while (bondedAtoms.hasNext()) {
					IAtom bondedAtom = bondedAtoms.next();
					if (!partContainer.contains(bondedAtom)
							&& !bondedAtom.getSymbol().equals("H"))
						partContainer.addAtom(bondedAtom);

				}
				partContainer.addBond(aBond);
			}

			partContainer.setProperty("Tipo", tipo);
			partContainer.setID(String.valueOf(ID));
			ID++;
			calcularCentroMasa(partContainer); // HAy que incluirle el H para el centro de masa
			centrosDescriptores.addAtomContainer(partContainer);

		}

		return centrosDescriptores;
	}

	private int numeroHeteroAtomos(IAtomContainer container){
		int n = 0;

		for(IAtom at : container.atoms()){
			if(!at.getSymbol().equals("H") && !at.getSymbol().equals("C"))
				n++;
		}

		return n;
	}






	/*
	 * Get
	 */

	public TreeMap<CDescrip, IAtomContainerSet> getCentrosDescriptores() {
		return cdescriptores;
	}
	
	public String getCaracteristicas(){
		String cad = "";
		//CDescrip[] desc = CDescrip.values();
		for(CDescrip tipo : CDescrip.values()){
			cad += tipo + "_" + cdescriptores.get(tipo).getAtomContainerCount() + "\n";
		}
		return cad;
	}

	public int size(){
		int temp = 0;
		for (Iterator<IAtomContainerSet> it= cdescriptores.values().iterator(); it.hasNext(); ) {
			IAtomContainerSet value = (IAtomContainerSet) it.next();
			temp += value.getAtomContainerCount();
		}
		return temp;
	}
	public IAtomContainer getMolecula() {
		return molecula;
	}

	public double[][] getDistanciasTopograficas(){
		return this.distanciasTopograficas;
	}
	public IAtomContainerSet getAnillos() {
		//return anillos;
		return cdescriptores.get(CDescrip.ANILLO);
	}

	public ArrayList<IAtom> getMetilo() {
		return metilo;
	}

	public IAtomContainerSet getMetiloSet() {
		//return metiloSet;
		return cdescriptores.get(CDescrip.METILO);
	}

	public ArrayList<IAtom> getMetileno() {
		return metileno;

	}

	public IAtomContainerSet getMetilenoSet() {
		//return metilenoSet;
		return cdescriptores.get(CDescrip.METILENO);
	}

	public ArrayList<IAtom> getMetino() {
		return metino;
	}

	public IAtomContainerSet getMetinoSet() {
		//return metinoSet;
		return cdescriptores.get(CDescrip.METINO);
	}

	public ArrayList<IAtom> getClusterOrden3() {

		return clusterOrden3;
	}

	public IAtomContainerSet getClusterOrden3Set() {

		//return clusterOrden3Set;
		return cdescriptores.get(CDescrip.CLUSTER3);
	}

	public ArrayList<IAtom> getClusterOrden4() {
		return clusterOrden4;
	}

	public IAtomContainerSet getClusterOrden4Set() {
		//return clusterOrden4Set;
		return cdescriptores.get(CDescrip.CLUSTER4);
	}

	public ArrayList<IAtom> getHeteroAtomo() {
		return heteroAtomo;
	}

	public IAtomContainerSet getHeteroAtomoSet() {
		//return heteroAtomoSet;
		return cdescriptores.get(CDescrip.HETEROATOMO);
	}

	public IAtomContainerSet getCentroDescriptor(CDescrip tipo){
		//		IAtomContainerSet temp = null;
		//		switch(tipo){
		//		case HETEROATOMO: temp = heteroAtomoSet; break;
		//		case CLUSTER3: temp = clusterOrden3Set; break;
		//		case CLUSTER4: temp = clusterOrden4Set; break;
		//		case METILO:  temp = metiloSet; break;
		//		case METINO: temp = metinoSet; break;
		//		case METILENO: temp = metilenoSet; break;
		//		case ANILLO: temp = anillos; break;
		//
		//		}
		//		return temp;

		return cdescriptores.get(tipo);
	}

	//	public ArrayList<List<List<IAtom>>> getCaminos(){
	//		
	//		this.caminos.trimToSize();
	//		return this.caminos;
	//	}
	public TreeMap <String, List<List<IAtom>>> getCaminos(){

		return this.caminos;
	}


	/*
	 * Set
	 */

	public void addCaminos(TreeMap <String, List<List<IAtom>>> c){

		this.caminos.putAll(c);
		//this.caminos.addAll(c);
		//this.caminos.trimToSize();

	}

	public void setDistanciaTopografica(double[][] dist){
		this.distanciasTopograficas = dist;
	}


	public void setHeteroAtomo(ArrayList<IAtom> heteroAtomo) {

		heteroAtomoSet = setVecinos(heteroAtomo, CDescrip.HETEROATOMO);
		//		for(IAtomContainer cont : heteroAtomoSet.atomContainers()){
		//			//calcularCentroMasa(cont);			
		//		}
		//CM(heteroAtomoSet);
		cdescriptores.put(CDescrip.HETEROATOMO, heteroAtomoSet);

		this.heteroAtomo = heteroAtomo;
	}

	public void setMolecula(IAtomContainer molecula) {
		this.molecula = molecula;
	}

	public void setAnillos(IAtomContainerSet anillos) {

		for(IAtomContainer at : anillos.atomContainers()){
			at.setProperty("Tipo", CDescrip.ANILLO);
			at.setID(String.valueOf(ID));
			ID++;
			calcularCentroMasa(at);			
		}
		//CM(anillos);
		cdescriptores.put(CDescrip.ANILLO, anillos);
	}

	public void setMetilo(ArrayList<IAtom> metilo) {
		metiloSet = setVecinos(metilo, CDescrip.METILO);
		//CM(metiloSet);
		cdescriptores.put(CDescrip.METILO, metiloSet);
		this.metilo = metilo;
	}

	public void setMetino(ArrayList<IAtom> metino) {
		metinoSet = setVecinos(metino, CDescrip.METINO);
		//CM(metinoSet);	
		cdescriptores.put(CDescrip.METINO, metinoSet);
		this.metino = metino;
	}

	public void setMetileno(ArrayList<IAtom> metileno) {
		metilenoSet = setVecinos(metileno, CDescrip.METILENO);
		//CM(metilenoSet);
		cdescriptores.put(CDescrip.METILENO, metilenoSet);
		this.metileno = metileno;
	}

	public void setClusterOrden3(ArrayList<IAtom> clusterOrden3) {
		clusterOrden3Set = setVecinos(clusterOrden3, CDescrip.CLUSTER3);
		for(@SuppressWarnings("unused") IAtomContainer cont : clusterOrden3Set.atomContainers()){
			//calcularCentroMasa(cont);			
		}
		//CM(clusterOrden3Set);
		cdescriptores.put(CDescrip.CLUSTER3, clusterOrden3Set);
		this.clusterOrden3 = clusterOrden3;
	}

	public void setClusterOrden4(ArrayList<IAtom> clusterOrden4) {
		clusterOrden4Set = setVecinos(clusterOrden4,CDescrip.CLUSTER4);
		//CM(clusterOrden4Set);
		cdescriptores.put(CDescrip.CLUSTER4, clusterOrden4Set);
		this.clusterOrden4 = clusterOrden4;
	}

	private void calcularCentroMasa(IAtomContainer container){
		Point3d cm = new Point3d();
		cm = container.getAtom(0).getPoint3d();
		String simb = container.getAtom(0).getSymbol();
		if(!simb.equals("F") &&  !simb.equals("Cl") && !simb.equals("Br") && !simb.equals("I")){
			cm = GeometryTools.get3DCentreOfMass(container);
		}
		//		else{
		//			cm = container.getAtom(0).getPoint3d();
		//			//Double mass = a.getExactMass();
		//			
		//			//cm = new Point3d(aux.x, aux.y/, aux.z/mass);
		//		}

		container.setProperty("CenterOfMass", cm);

	}

	//	private void CM(IAtomContainerSet at){
	//		for(IAtomContainer cont : at.atomContainers()){
	//			System.out.println(cont.getProperty("Tipo").toString());
	//			System.out.println(cont.getProperty("CenterOfMass").toString());			
	//		}
	//	}

}
