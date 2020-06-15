package domain;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Point3d;

import domain.ContenedorCD.CDescrip;

import org.omegahat.Environment.DataStructures.list;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.UnsupportedChemObjectException;
import org.openscience.cdk.geometry.GeometryTools;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topographical.ElectroTopographicDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topographical.LipoTopographicDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topographical.RefractedTopographicDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topological.ElectroTopologicalDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topological.LipoTopologicDescriptor;
import org.openscience.cdk.qsar.descriptors.hybrid.topological.RefractedTopologicalDescriptor;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.RingPartitioner;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.similarity.Cosine;
import org.openscience.cdk.similarity.Jaccard;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.tools.LoggingTool;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

//import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader.Array;

import util.Data;
import util.DescriptorCenterTable;
import util.ExcelStadistics;
import util.GruposFuncionalesException;
import util.MatrixCompared;
import util.MoleculaException;
import util.Normalization;

@SuppressWarnings({ "unused" })
public class Molecule {

	private final class MDLReaderExtension extends MDLReader {
		private MDLReaderExtension(Reader in) {
			super(in);
		}

		public int read(char[] arg0, int arg1, int arg2) throws IOException {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public void close() throws IOException {
			// TODO Auto-generated method stub
			
		}
	}


	// private static Molecule contCentroDescriptor = new Molecule();
	IAtomContainer graph;
	int[][] matrizTopologica;
	private String direccion;
	private String nombre;
	private String actividad;
	private ReducedGraph reducedGraph;
	private ArrayList<IAtom> listAtomSelected;
	private double[][] distancia;
	SimilaryFunctions similaryFunction;
	private NumberFormat format = NumberFormat.getInstance();
	private Normalization normalization;

	public Molecule() {
		similaryFunction = new SimilaryFunctions();
		reducedGraph = new ReducedGraph();
		listAtomSelected = new ArrayList<IAtom>();
		format.setMaximumFractionDigits(5);
		normalization = null;
	}

	public void setMolecula(String file) throws FileNotFoundException, CDKException {
		MDLReader reader = new MDLReader(new FileReader(new File(file)));
		ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
		List<IAtomContainer> containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
		IAtomContainer molecule = containersList.get(0);
		this.graph = molecule;
	}
	
	public void setNormalization(Normalization normalization){
		this.normalization = normalization;
	}

	public void setMolecula(File file) throws FileNotFoundException, CDKException {
		try {
			direccion = file.getAbsolutePath();
			actividad = "";
			MDLReader reader = new MDLReader(new FileReader(file));			
			ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			List<IAtomContainer> containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
			IAtomContainer molecula = containersList.get(0);
			molecula.setID(file.getName());
			for(int i=0; i< molecula.getAtomCount(); i++)
				molecula.getAtom(i).setID(String.valueOf(i+1));
			nombre = file.getName();
			graph = molecula;
		} catch (Exception e) {
			throw new MoleculaException("Molecula no valida");
		}
	}

	public ArrayList<IAtom> getAtomsSelected() {
		return listAtomSelected;
	}
	
	public Normalization getNormalization(){
		return normalization;
	}

	public void setAtomsSelected(IAtom atom) {
		listAtomSelected.add(atom);
	}

	public void clearAtomSelected() {
		listAtomSelected.clear();
	}

	public void setActividad(String valor) {
		actividad = valor;
	}

	public String getActividad() {
		return actividad;
	}

	public String getNombre() {
		return nombre;
	}

	public String getDireccion() {
		return direccion;
	}

	public void ReducedGraph() throws FileNotFoundException, GruposFuncionalesException, UnsupportedChemObjectException {
		reducedGraph = buscarAnillos(reducedGraph);
		reducedGraph = buscarCluster(reducedGraph);
		reducedGraph = buscarGrupoFuncional3(reducedGraph);
		reducedGraph = buscarHetroatomos(reducedGraph);
		
		calcularDistTopograficas(reducedGraph);
	}	
	
	public ReducedGraph getReducedGraph() throws FileNotFoundException, GruposFuncionalesException, UnsupportedChemObjectException {
		reducedGraph = buscarAnillos(reducedGraph);
		reducedGraph = buscarCluster(reducedGraph);
		reducedGraph = buscarGrupoFuncional3(reducedGraph);
		//reducedGraph = buscarHetroatomos(reducedGraph);
		@SuppressWarnings("rawtypes")
		Iterator cd = reducedGraph.getAllCentrosDescriptores().iterator();
		while (cd.hasNext()) {
			CentroDescriptor e = (CentroDescriptor) cd.next();
			e.AllIndiceTotalTopographical();
		}
		return reducedGraph;
	}
	
	public ContenedorCD buscarGrupoFuncional(ContenedorCD centrosDescriptores) throws CDKException, FileNotFoundException {
		ArrayList<IAtom> heteroAtomos = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetilo = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetileno = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetino = new ArrayList<IAtom>();
		ArrayList<IAtom> clusterOrden3 = new ArrayList<IAtom>();
		ArrayList<IAtom> clusterOrden4 = new ArrayList<IAtom>();
		for (IAtom atomo : this.graph.atoms()) {
			if (!atomo.getSymbol().equals("H") && !Boolean.parseBoolean(atomo.getProperty("CD").toString())) {
				List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomo);
				int hidrogenos = numeroHidrogenos(atomosConectados);
				// Busca los diferentes centros descriptores excepto los anillos
				IBond.Order c = graph.getMaximumBondOrder(atomo);
				boolean metilo = false, metileno = false, metino = false, cluster3 = false, cluster4 = false;
				if (hidrogenos == 3 && atomo.getSymbol().equals("C"))
					metilo = true;
				else if (hidrogenos == 2 && atomo.getSymbol().equals("C") && c.equals(IBond.Order.DOUBLE))
					metileno = true;
				else if (hidrogenos == 1 && atomo.getSymbol().equals("C") && c.equals(IBond.Order.TRIPLE))
					metino = true;
				else {
					int numeroAtomosConectados = graph.getConnectedAtomsCount(atomo);
					if ((hidrogenos == 1 && numeroAtomosConectados == 4) || (hidrogenos == 0 && numeroAtomosConectados == 3))
						cluster3 = true;
					else if (hidrogenos == 0 && numeroAtomosConectados == 4)
						cluster4 = true;
				}
				if (metilo)
					grupoMetilo.add(atomo); // Guardar en cola prioridad
				else if (metileno)
					grupoMetileno.add(atomo); // Guardar en cola prioridad
				else if (metino)
					grupoMetino.add(atomo); // Guardar en cola prioridad
				else if (cluster3)
					clusterOrden3.add(atomo);
				else if (cluster4)
					clusterOrden4.add(atomo);
				if (!atomo.getSymbol().equals("C")) // Guardar en cola prioridad
					heteroAtomos.add(atomo);
			}
		}
		// ContenedorCD centrosDescriptores = new ContenedorCD(this.molecula);
		centrosDescriptores.setClusterOrden3(clusterOrden3);
		centrosDescriptores.setClusterOrden4(clusterOrden4);
		centrosDescriptores.setHeteroAtomo(heteroAtomos);
		centrosDescriptores.setMetileno(grupoMetileno);
		centrosDescriptores.setMetilo(grupoMetilo);
		centrosDescriptores.setMetino(grupoMetino);
		return centrosDescriptores;
	}

	public void selectAtoms(CentroDescriptor cd) {
		if (cd.getTipo() == TipoCentroDescriptor.ANILLO || cd.getTipo() == TipoCentroDescriptor.CLUSTER3 || cd.getTipo() == TipoCentroDescriptor.CLUSTER4) {
			for (int i = 0; i < cd.getCantidadAtomos(); i++)
				for (int j = 0; j < getCantAtmos(); j++)
					if (getGraph().getAtom(j).getID().equals(cd.getIdAtomo(i))) {
						getGraph().getAtom(j).setFlag(0, true);
						break;
					}
		} else {
			for (int j = 0; j < getCantAtmos(); j++)
				if (getGraph().getAtom(j).getID().equals(cd.getIdAtomo(0))) {
					getGraph().getAtom(j).setFlag(0, true);
					break;
				}
		}
	}

	private int numeroHidrogenos(List<IAtom> atomosConectados) {
		int cHidrogeno = 0;
		for (Iterator<IAtom> atomo = atomosConectados.iterator(); atomo.hasNext();) {
			IAtom iAtom = (IAtom) atomo.next();
			if (iAtom.getSymbol().equals("H"))
				cHidrogeno++;
		}
		return cHidrogeno;
	}

	public ContenedorCD buscarGrupoFuncional2(ContenedorCD centrosDescriptores) throws CDKException, FileNotFoundException {
		ArrayList<IAtom> heteroAtomos = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetilo = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetileno = new ArrayList<IAtom>();
		ArrayList<IAtom> grupoMetino = new ArrayList<IAtom>();
		ArrayList<IAtom> clusterOrden3 = new ArrayList<IAtom>();
		ArrayList<IAtom> clusterOrden4 = new ArrayList<IAtom>();
		for (IAtom atomo : this.graph.atoms())
			if (!atomo.getSymbol().equals("H") && !Boolean.parseBoolean(atomo.getProperty("CD").toString())) {
				CDescrip tipoAtomoCentral = Tipo(atomo);
				// String n = atomo.getID();
				if (tipoAtomoCentral == CDescrip.CLUSTER4) {
					Marcar(atomo);
					// System.out.println(n);
					clusterOrden4.add(atomo);
				} else if (tipoAtomoCentral == CDescrip.CLUSTER3) {
					Marcar(atomo);
					// System.out.println(n);
					clusterOrden3.add(atomo);
				} else if (tipoAtomoCentral.compareTo(CDescrip.METINO) >= 0 && tipoAtomoCentral.compareTo(CDescrip.HETEROATOMO) <= 0) {
					if (!ConectadoCluster(atomo)) {
						atomo.setProperty("CD", true);
						if (tipoAtomoCentral.compareTo(CDescrip.METINO) == 0) {
							grupoMetino.add(atomo);
							// System.out.println(n);
						} else if (tipoAtomoCentral.compareTo(CDescrip.METILO) == 0) {
							grupoMetilo.add(atomo);
							// System.out.println(n);
						} else if (tipoAtomoCentral.compareTo(CDescrip.METILENO) == 0) {
							// System.out.println(n);
							grupoMetileno.add(atomo);
						} else if (tipoAtomoCentral.compareTo(CDescrip.HETEROATOMO) == 0) {
							// System.out.println(n);
							heteroAtomos.add(atomo);
						}
					}
				}
			}
		// ContenedorCD centrosDescriptores = new ContenedorCD(this.molecula);
		centrosDescriptores.setClusterOrden3(clusterOrden3);
		centrosDescriptores.setClusterOrden4(clusterOrden4);
		centrosDescriptores.setHeteroAtomo(heteroAtomos);
		centrosDescriptores.setMetileno(grupoMetileno);
		centrosDescriptores.setMetilo(grupoMetilo);
		centrosDescriptores.setMetino(grupoMetino);
		return centrosDescriptores;
	}

	private CDescrip Tipo(IAtom atomoCentral) {
		CDescrip tipo = CDescrip.ANILLO;
		List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomoCentral);
		int hidrogenos = numeroHidrogenos(atomosConectados);
		// Busca los diferentes centros descriptores excepto los anillos
		IBond.Order c = graph.getMaximumBondOrder(atomoCentral);
		// boolean metilo = false, metileno = false, metino = false, cluster3 =
		// false, cluster4 = false;
		if (hidrogenos == 3 && atomoCentral.getSymbol().equals("C"))
			tipo = CDescrip.METILO;
		else if (hidrogenos == 2 && atomoCentral.getSymbol().equals("C") && c.equals(IBond.Order.DOUBLE))
			tipo = CDescrip.METILENO;
		else if (hidrogenos == 1 && atomoCentral.getSymbol().equals("C") && c.equals(IBond.Order.TRIPLE))
			tipo = CDescrip.METINO;
		else if (!atomoCentral.getSymbol().equals("C"))
			tipo = CDescrip.HETEROATOMO;
		else {
			int numeroAtomosConectados = graph.getConnectedAtomsCount(atomoCentral);
			if ((hidrogenos == 1 && numeroAtomosConectados == 4) || (hidrogenos == 0 && numeroAtomosConectados == 3))
				tipo = CDescrip.CLUSTER3;
			else if (hidrogenos == 0 && numeroAtomosConectados == 4)
				tipo = CDescrip.CLUSTER4;
		}
		return tipo;
	}

	private void Marcar(IAtom atomoCentral) {
		List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomoCentral);
		for (Iterator<IAtom> atomo = atomosConectados.iterator(); atomo.hasNext();) {
			IAtom iAtom = (IAtom) atomo.next();
			if (!iAtom.getSymbol().equals("H")) {
				// iAtom.setProperty("CD", true);
				iAtom.setFlag(11, true);
				if(!iAtom.getSymbol().equals("C"))
					iAtom.setFlag(12, true);
			}
		}
		// atomoCentral.setProperty("CD", true);
		atomoCentral.setFlag(11, true);
	}

	private boolean ConectadoCluster(IAtom atomoCentral) {
		boolean result = false;
		List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomoCentral);
		// CDescrip mayor = CDescrip.HETEROATOMO;
		TipoCentroDescriptor mayor = TipoCentroDescriptor.HETEROATOMO;
		for (Iterator<IAtom> atomo = atomosConectados.iterator(); atomo.hasNext();) {
			IAtom iAtom = (IAtom) atomo.next();
			if (!iAtom.getSymbol().equals("H")) {
				// CDescrip aux = Tipo(iAtom);
				TipoCentroDescriptor aux = Tipo2(iAtom);
				// if(mayor.compareTo(aux) < 0)
				if (mayor.compareTo(aux) < 0)
					mayor = aux;
			}
		}
		// if(mayor.compareTo(CDescrip.CLUSTER3) == 0 ||
		// mayor.compareTo(CDescrip.CLUSTER4)== 0)
		if (mayor == TipoCentroDescriptor.CLUSTER3 || mayor == TipoCentroDescriptor.CLUSTER4)
			result = true;
		return result;
	}

	public ReducedGraph buscarGrupoFuncional3(ReducedGraph reducedGraph) throws FileNotFoundException, UnsupportedChemObjectException, GruposFuncionalesException {
		try {
			for (IAtom atomo : graph.atoms()) {
				checkValidAtom(atomo);
				if (!atomo.getSymbol().equals("H") && !atomo.getFlag(11)) {
					TipoCentroDescriptor tipoAtomoCentral = Tipo2(atomo);
					IAtomContainer aux = null;
					if (!tipoAtomoCentral.equals(TipoCentroDescriptor.NINGUNO)) {
						switch (tipoAtomoCentral) {
							case METILENO:
							case METILO:
							case METINO:
							case CARBONOCONECTADO:
							case HETEROATOMO:
								if (!ConectadoCluster2(atomo, tipoAtomoCentral)) {
									atomo.setFlag(11, true);
									aux = setVecinos(atomo, tipoAtomoCentral);
								}
								break;
						default:
							break;
						}
						if (aux != null)
							reducedGraph.addCentroDescriptor(new CentroDescriptor(aux, tipoAtomoCentral));
					}
				}
			}
		} catch (Exception e) {
			throw new GruposFuncionalesException("Otros");
		}
		return reducedGraph;
	}
	
	public ReducedGraph buscarHetroatomos(ReducedGraph reducedGraph) throws FileNotFoundException, UnsupportedChemObjectException, GruposFuncionalesException {
		try {
			for (IAtom atomo : graph.atoms()) {
				checkValidAtom(atomo);
				if ((!atomo.getSymbol().equals("H") && !atomo.getSymbol().equals("C")) && (!atomo.getFlag(11) || atomo.getFlag(12))) {
					TipoCentroDescriptor tipoAtomoCentral = Tipo2(atomo);
					IAtomContainer aux = null;
					if (!tipoAtomoCentral.equals(TipoCentroDescriptor.NINGUNO)) {
						switch (tipoAtomoCentral) {
							case HETEROATOMO:
								if (!ConectadoCluster2(atomo, tipoAtomoCentral)) {
									atomo.setFlag(11, true);
									aux = setVecinos(atomo, tipoAtomoCentral);
								}
								break;
						default:
							break;
						}
						if (aux != null)
							reducedGraph.addCentroDescriptor(new CentroDescriptor(aux, tipoAtomoCentral));
					}
				}
			}
		} catch (Exception e) {
			throw new GruposFuncionalesException("Otros");
		}
		return reducedGraph;
	}

	private boolean ConectadoCluster2(IAtom atomoCentral, TipoCentroDescriptor tipoAtomoCentral) {
		// boolean result = false;
		List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomoCentral);
		// CDescrip mayor = CDescrip.HETEROATOMO;
		TipoCentroDescriptor mayor = tipoAtomoCentral;
		for (Iterator<IAtom> atomo = atomosConectados.iterator(); atomo.hasNext();) {
			IAtom iAtom = (IAtom) atomo.next();
			if (!iAtom.getSymbol().equals("H")) {
				// CDescrip aux = Tipo(iAtom);
				TipoCentroDescriptor aux = Tipo2(iAtom);

				if ((aux == TipoCentroDescriptor.CLUSTER3 || aux == TipoCentroDescriptor.CLUSTER4) && !iAtom.getFlag(11))
					mayor = aux;
			}
		}
		// if(mayor.compareTo(CDescrip.CLUSTER3) == 0 ||
		// mayor.compareTo(CDescrip.CLUSTER4)== 0)
		// if(mayor == TipoCentroDescriptor.CLUSTER3 || mayor ==
		// TipoCentroDescriptor.CLUSTER4)
		// result = true;
		// return result;
		return (mayor.equals(TipoCentroDescriptor.CLUSTER3) || mayor.equals(TipoCentroDescriptor.CLUSTER4));
	}

	private void checkValidAtom(IAtom atomo) throws UnsupportedChemObjectException {
		// Cambiar equals por contains
		if (!atomo.getSymbol().equals("C") && !atomo.getSymbol().equals("H") && !atomo.getSymbol().equals("O") && !atomo.getSymbol().equals("S") && !atomo.getSymbol().equals("N") && !atomo.getSymbol().equals("P") && !atomo.getSymbol().equals("F") && !atomo.getSymbol().equals("Cl") && !atomo.getSymbol().equals("Br") && !atomo.getSymbol().equals("I") /*
																																																																																								 * && ! atomo . getSymbol ( ) . equals ( "As" )																																																																																					 */)
			throw new UnsupportedChemObjectException("Elemento no valido");
	}

	private TipoCentroDescriptor Tipo2(IAtom atomoCentral) {
		TipoCentroDescriptor tipo = TipoCentroDescriptor.NINGUNO;
		List<IAtom> atomosConectados = graph.getConnectedAtomsList(atomoCentral);
		int hidrogenos = numeroHidrogenos(atomosConectados);
		int numeroAtomosConectados = graph.getConnectedAtomsCount(atomoCentral);
		// Busca los diferentes centros descriptores excepto los anillos
		IBond.Order c = graph.getMaximumBondOrder(atomoCentral);
		// boolean metilo = false, metileno = false, metino = false, cluster3 =
		// false, cluster4 = false;
		if (hidrogenos == 3 && atomoCentral.getSymbol().equals("C"))
			tipo = TipoCentroDescriptor.METILO;
		else if (hidrogenos == 2 && atomoCentral.getSymbol().equals("C") && c.equals(IBond.Order.DOUBLE))
			tipo = TipoCentroDescriptor.METILENO;
		else if (hidrogenos == 1 && atomoCentral.getSymbol().equals("C") && c.equals(IBond.Order.TRIPLE))
			tipo = TipoCentroDescriptor.METINO;
		else if ((hidrogenos == 1 && numeroAtomosConectados == 4) || (hidrogenos == 0 && numeroAtomosConectados == 3))
				tipo = TipoCentroDescriptor.CLUSTER3;
		else if (hidrogenos == 0 && numeroAtomosConectados == 4)
				tipo = TipoCentroDescriptor.CLUSTER4;
		else if (!atomoCentral.getSymbol().equals("C"))
			tipo = TipoCentroDescriptor.HETEROATOMO;
		else if ((hidrogenos == 2 && numeroAtomosConectados == 4)|| (hidrogenos == 1 && numeroAtomosConectados == 3) || (numeroAtomosConectados == 2 && c.equals(IBond.Order.TRIPLE)))
		    tipo = TipoCentroDescriptor.CARBONOCONECTADO;

		return tipo;
	}

	private IAtomContainer setVecinos(IAtom atomocentral, TipoCentroDescriptor tipo) {
		ArrayList<IBond> enlaces = new ArrayList<IBond>(graph.getConnectedBondsList(atomocentral));
		IAtomContainer partContainer = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
		partContainer.addAtom(atomocentral);
		for (IBond aBond : enlaces) {
			Iterator<IAtom> bondedAtoms = aBond.atoms().iterator();
			while (bondedAtoms.hasNext()) {
				IAtom bondedAtom = bondedAtoms.next();
				if (tipo.equals(TipoCentroDescriptor.HETEROATOMO)){
					if (bondedAtom.getSymbol().equals("H"))
					    if (!partContainer.contains(bondedAtom))
					        partContainer.addAtom(bondedAtom);
				}else {
				   if (!partContainer.contains(bondedAtom))
						   partContainer.addAtom(bondedAtom);
				}
			}
			partContainer.addBond(aBond);
		}
		// partContainer.setProperty("Tipo", tipo);
		// partContainer.setID(String.valueOf(ID));
		// ID++;
		// calcularCentroMasa(partContainer); // HAy que incluirle el H para el
		// centro de masa
		// centrosDescriptores.addAtomContainer(partContainer);
		return partContainer;
	}

	public ReducedGraph buscarAnillos(ReducedGraph reducedGraph) throws FileNotFoundException, GruposFuncionalesException {
		try {
			//Buscar anillos menores que 10.
			@SuppressWarnings("unchecked")
			List<IRingSet> ciclosAux = new SSSRFinder(graph).findEquivalenceClasses();
			for (int f = 0; f < ciclosAux.size(); f++) {
				IAtomContainer ac = RingPartitioner.convertToAtomContainer(ciclosAux.get(f));
				for (IAtom atom : ac.atoms())
					if (!atom.getSymbol().equals("H"))
						atom.setFlag(11, true);
				CentroDescriptor cd = new CentroDescriptor(ac, TipoCentroDescriptor.ANILLO);
				reducedGraph.addCentroDescriptor(cd);
			}
			
			//Buscar anillos mayores que 10 y menores que 13.
			ArrayList<IAtomContainer> listAC = new ArrayList<IAtomContainer>();
			for (int i = 0; i < ciclosAux.size()-1; i++)
				for (int j = i+1; j < ciclosAux.size(); j++) {
					IAtomContainer temp = anillosCompuestos(ciclosAux.get(i), ciclosAux.get(j));
					if(temp != null) {
						if (listAC.isEmpty())
							listAC.add(temp);
						else {
							boolean flag = false;
							for(int l = 0; l < listAC.size(); l++)
								if(ringEquals(temp, listAC.get(l))){
									flag = true;
									break;
								}
							if(!flag)
								listAC.add(temp);
						}
					}
				}
			for (IAtomContainer iAtomContainer : listAC) {
				CentroDescriptor cd = new CentroDescriptor(iAtomContainer, TipoCentroDescriptor.ANILLO);
			    reducedGraph.addCentroDescriptor(cd);
			}	
			
			//Buscar anillos mayores que 14 y menores de 17.
			ArrayList<IAtomContainer> listAC1 = new ArrayList<IAtomContainer>();
			if (listAC.size() > 1) {
				for (int i = 0; i < listAC.size() - 1; i++)
					for(int j = i + 1; j < listAC.size(); j++) {
						IAtomContainer temp = anillosCompuestosIAtomContainer(listAC.get(i), listAC.get(j));
						if(temp != null)
							if (listAC1.isEmpty())
								listAC1.add(temp);
							else {
								boolean flag = false;
								for (int l = 0; l < listAC1.size(); l++)
									if (ringEquals(temp, listAC1.get(l))){
										flag = true;
										break;
									}
								if(!flag)
									listAC1.add(temp);
							}
					}
				for (IAtomContainer iAtomContainer : listAC1) {
					CentroDescriptor cd = new CentroDescriptor(iAtomContainer, TipoCentroDescriptor.ANILLO);
				    reducedGraph.addCentroDescriptor(cd);
				}
			}	
			
			//Buscar anillos mayores que 18.
			ArrayList<IAtomContainer> listAC2 = new ArrayList<IAtomContainer>();
			if (listAC1.size() > 1) {
				int mayor = listAC1.get(0).getAtomCount();
				for (int i = 1; i < listAC1.size() - 1; i++)
					if(mayor < listAC1.get(i).getAtomCount())
						mayor = listAC1.get(i).getAtomCount();
				for (int i = 0; i < listAC1.size() - 1; i++)
					for(int j = i + 1; j < listAC1.size(); j++) {
						IAtomContainer temp = anillosCompuestosIAtomContainer(listAC1.get(i), listAC1.get(j));
						if(temp != null)
							if (listAC2.isEmpty())
								listAC2.add(temp);
							else {
								boolean flag = false;
								for (int l = 0; l < listAC2.size(); l++)
									if (ringEquals(temp, listAC2.get(l))){
										flag = true;
										break;
									}
								if(!flag && temp.getAtomCount() > mayor)
									listAC2.add(temp);
							}
					}
				for (IAtomContainer iAtomContainer : listAC2) {
					CentroDescriptor cd = new CentroDescriptor(iAtomContainer, TipoCentroDescriptor.ANILLO);
				    reducedGraph.addCentroDescriptor(cd);
				}
			}	
		} catch (Exception e) {
			throw new GruposFuncionalesException("Anillo");
		}
		return reducedGraph;
	}
	
	public boolean ringEquals(IAtomContainer ring1, IAtomContainer ring2){
		boolean flag = false; int cont =0;
		if(ring1.getAtomCount()==ring2.getAtomCount()){
			for (IAtom atom : ring1.atoms())
				for (IAtom atom1 : ring2.atoms()) 
					if(atom.getSymbol().equals(atom1.getSymbol()) && atom.getID().equals(atom1.getID())) {
						cont ++;
					} 
			if(cont == ring1.getAtomCount())
				flag = true;
		}
		return flag;
	}
	
    public ReducedGraph buscarCluster(ReducedGraph reducedGraph) throws FileNotFoundException, UnsupportedChemObjectException, GruposFuncionalesException {
    	
    	for (IAtom atomo : graph.atoms()) {
    		checkValidAtom(atomo);
			if (!atomo.getSymbol().equals("H")) {
				TipoCentroDescriptor tipoAtomoCentral = Tipo2(atomo);
				IAtomContainer aux = null;
				if (tipoAtomoCentral.equals(TipoCentroDescriptor.CLUSTER3)||tipoAtomoCentral.equals(TipoCentroDescriptor.CLUSTER4))
					if(!cicleAtomContainer(atomo)){
						Marcar(atomo);
						aux = setVecinos(atomo, tipoAtomoCentral);
					}
				if (aux != null)
					reducedGraph.addCentroDescriptor(new CentroDescriptor(aux, tipoAtomoCentral));
			}
    	}
    	
    	return reducedGraph;
    }
     
	public boolean cicleAtomContainer(IAtom atomo){
		boolean flag = false;
		
		@SuppressWarnings("unchecked")
		Iterator<CentroDescriptor> listCentro = reducedGraph.getCentDescripPorTipo(TipoCentroDescriptor.ANILLO).iterator();
		while (listCentro.hasNext()){
			CentroDescriptor cd = listCentro.next();
			for(int i=0; i< cd.getCantidadAtomos(); i++){
			   if(cd.getSimbolAtomo(i).equals(atomo.getSymbol()) && cd.getIdAtomo(i).equals(atomo.getID()))	
				   return flag = true;
			}
				
		}
		
		return flag;
	}
	
	public boolean equalsIAtomContainer (IAtomContainer cd1, IAtomContainer cd2) {
		boolean flag = false;
		int cont = 0;
		
		for (IAtom atom : cd1.atoms()){
			String b1 = atom.getSymbol() + atom.getID();
			for (IAtom atom1 : cd2.atoms()) {
				String b2 = atom1.getSymbol() + atom1.getID();
				if(b1.equals(b2))
					cont ++;
				//break;
			}
		}
		
		if(cont == cd2.getAtomCount()||cont == cd1.getAtomCount())
			flag = true;
	
		return flag;
	}
	
	public IAtomContainer anillosCompuestos(IRingSet cd1, IRingSet cd2){
		int cont = 0;
		IAtomContainer ac = null;
		for (IAtom atom : cd1.getAtomContainer(0).atoms()) 
			for (IAtom atom1 : cd2.getAtomContainer(0).atoms()) 
				if(atom.getSymbol().equals(atom1.getSymbol()) && atom.getID().equals(atom1.getID())) {
					cont ++;
				}
			
		if(cont == 2 || cont == 3 ) {
			ac = RingPartitioner.convertToAtomContainer(cd1);
			IAtomContainer ac1 = RingPartitioner.convertToAtomContainer(cd2);
			ac.add(ac1);
		}
		return ac;
	}
	
	public IAtomContainer anillosCompuestosIAtomContainer(IAtomContainer cd1, IAtomContainer cd2){
		int cont = 0;
		ArrayList<IAtom> listAtom = new ArrayList<IAtom>();
		for (IAtom atom : cd1.atoms())
			for (IAtom atom1 : cd2.atoms()) 
				if(atom.getSymbol().equals(atom1.getSymbol()) && atom.getID().equals(atom1.getID())) {
					cont ++;
					listAtom.add(atom1);
				} 
		
		if(cont >= 2) {
			IAtomContainer ac = new AtomContainer();
			ac.add(cd1);
			for(IAtom atom : cd2.atoms())
				for (IAtom atom1 : listAtom)
					if(!(atom.getSymbol().equals(atom1.getSymbol()) && atom.getID().equals(atom1.getID())))
				          ac.addAtom(atom);	
			return ac;
		}
		return null;
	}

	public int getAllCaminos(ContenedorCD centros) {
		int suma = 0;
		CDescrip[] desc = CDescrip.values();
		for (int i = 0; i < desc.length; i++) 
			for (int j = i; j < desc.length; j++) 
				suma += getCaminos(desc[i], desc[j], centros);
		return suma;
	}

	public int getCaminos(CDescrip tipo1, CDescrip tipo2, ContenedorCD centros) {
		IAtomContainerSet container1 = centros.getCentroDescriptor(tipo1);
		IAtomContainerSet container2 = centros.getCentroDescriptor(tipo2);
		int cant1 = container1.getAtomContainerCount();
		int cant2 = container2.getAtomContainerCount();
		for (int i = 0; i < cant1; i++) {
			IAtomContainer mol = container1.getAtomContainer(i);
			for (int j = 0; j < cant2; j++) {
				IAtomContainer mol2 = container2.getAtomContainer(j);
				TreeMap<String, List<List<IAtom>>> temp = BuscarCaminos(mol, mol2);
				centros.addCaminos(temp);
			}
		}
		return cant1 * cant2;
	}

	private TreeMap<String, List<List<IAtom>>> BuscarCaminos(IAtomContainer cont1, IAtomContainer cont2) {
		TreeMap<String, List<List<IAtom>>> map = new TreeMap<String, List<List<IAtom>>>();
		// ArrayList<List<List<IAtom>>> caminos = new
		// ArrayList<List<List<IAtom>>>();
		// int i, j;
		// i = j = 0;
		for (IAtom a : cont1.atoms()) {
			// j++;
			// i = 0;
			for (IAtom b : cont2.atoms()) {
				List<List<IAtom>> caminoSimple = PathTools.getAllPathsCalede(this.graph, a, b);
				String key = this.graph.getID() + " " + a.getSymbol() + a.getID() + "to" + b.getSymbol() + b.getID();
				map.put(key, caminoSimple);
			}
		}
		return map;
	}

	public void calcularDistTopograficas(ReducedGraph reducedGraph) {
		double dist = 0.0;
		int n = reducedGraph.getNumeroCentrosDescriptores();
		Object[] cd = reducedGraph.getAllCentrosDescriptores().toArray();
		distancia = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				CentroDescriptor cd1 = (CentroDescriptor) cd[i];
				CentroDescriptor cd2 = (CentroDescriptor) cd[j];
				boolean tt = includeCentroDescriptor(cd1, cd2);
				if (i == j || tt)
					dist = 0.0;
				else{
					dist = cd1.getCentroMasa().distance(cd2.getCentroMasa());
				}
				// reducedGraph.addCaminoTopografico(cd1, cd2, dist);
				distancia[i][j] = dist;
			}
		}
	}
	
	public void calcularDistTopograficas(Molecule molecule) {
		double dist = 0.0;
		Object[] atoms = molecule.getAtoms().toArray();
		distancia = new double[atoms.length][atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			for (int j = 0; j < atoms.length; j++) {
				String atom1 = ((IAtom)atoms[i]).getID() + ((IAtom)atoms[i]).getSymbol();
				String atom2 = ((IAtom)atoms[j]).getID() + ((IAtom)atoms[j]).getSymbol();
				boolean tt = atom1==atom2;
				if (i == j || tt)
					dist = 0.0;
				else{
					dist = ((IAtom)atoms[i]).getPoint3d().distance(((IAtom)atoms[j]).getPoint3d());
				}
				// reducedGraph.addCaminoTopografico(cd1, cd2, dist);
				//System.out.println(i + "-"+ j + "-" +dist);
				distancia[i][j] = dist;
			}
		}
	}
	
	private boolean includeCentroDescriptor (CentroDescriptor cd1, CentroDescriptor cd2) {
		boolean flag = false;
		int cantAtomos, cont = 0;

		if(cd1.getCantidadAtomos() < cd2.getCantidadAtomos())
			cantAtomos = cd1.getCantidadAtomos();
		else
			cantAtomos = cd2.getCantidadAtomos();
		
		for(int i=0; i < cd1.getCantidadAtomos(); i++){
			String f1= cd1.getSimbolAtomo(i) + cd1.getIdAtomo(i);
			for(int j=0; j < cd2.getCantidadAtomos(); j++) { 
				String f2= cd2.getSimbolAtomo(j) + cd2.getIdAtomo(j);
				if(f1.equals(f2)){
					cont ++;
					break;
				}	
			}	
		}
		
		if(cont == cantAtomos)
			flag = true;
		
		return flag;
	}

	double[][] dist;

	public void calcularTodasDistanciasTopograficas(ContenedorCD centros) {
		int numCD = centros.size();
		dist = new double[numCD][numCD];
		double fillValue = 0;
		for (int k = 0; k < numCD; k++)
			Arrays.fill(dist[k], fillValue);
		CDescrip[] desc = CDescrip.values();
		for (int i = 0; i < desc.length; i++) {
			for (int j = i; j < desc.length; j++) {
				if (i != j)
					getDistanciasTopograficasCDDiferentes(desc[i], desc[j], centros, dist);
				else
					getDistanciasTopograficasCDIguales(desc[i], centros, dist);
			}
		}
		//int a = dist.length;
		centros.setDistanciaTopografica(dist);
	}

	private void getDistanciasTopograficasCDDiferentes(CDescrip tipo1, CDescrip tipo2, ContenedorCD centros, double[][] dist) {
		IAtomContainerSet container1 = centros.getCentroDescriptor(tipo1);
		// for(IAtomContainer at : container1.atomContainers())
		// System.out.println(at.getID());
		IAtomContainerSet container2 = centros.getCentroDescriptor(tipo2);
		int cant1 = container1.getAtomContainerCount();
		int cant2 = container2.getAtomContainerCount();
		for (int i = 0; i < cant1; i++) {
			IAtomContainer mol = container1.getAtomContainer(i);
			for (int j = 0; j < cant2; j++) {
				IAtomContainer mol2 = container2.getAtomContainer(j);
				double value = DistanciaTopografica(mol, mol2);
				dist[Integer.valueOf(mol.getID())][Integer.valueOf(mol2.getID())] = value;
				dist[Integer.valueOf(mol2.getID())][Integer.valueOf(mol.getID())] = value;
				// System.out.println(value);
			}
		}
	}

	private void getDistanciasTopograficasCDIguales(CDescrip tipo, ContenedorCD centros, double[][] dist) {
		IAtomContainerSet container = centros.getCentroDescriptor(tipo);
		// for(IAtomContainer at : container1.atomContainers())
		// System.out.println(at.getID());
		int cant = container.getAtomContainerCount();
		for (int i = 0; i < cant; i++) {
			IAtomContainer mol = container.getAtomContainer(i);
			for (int j = i + 1; j < cant; j++) {
				IAtomContainer mol2 = container.getAtomContainer(j);
				double value = DistanciaTopografica(mol, mol2);
				dist[Integer.valueOf(mol.getID())][Integer.valueOf(mol2.getID())] = value;
				dist[Integer.valueOf(mol2.getID())][Integer.valueOf(mol.getID())] = value;
				// System.out.println(value);
			}
		}
	}

	private double DistanciaTopografica(IAtomContainer mol, IAtomContainer mol2) {
		Point3d p = (Point3d) mol.getProperty("CenterOfMass");
		Point3d p2 = (Point3d) mol2.getProperty("CenterOfMass");
		return p.distance(p2);
	}

	public void calcularIndicesTopologicos() throws CDKException, CloneNotSupportedException {
		DescriptorValue value;
		IAtomContainer aux = (IAtomContainer) this.graph.clone();
		IHybridDescriptor refract = new RefractedTopologicalDescriptor();
		IHybridDescriptor elect = new ElectroTopologicalDescriptor();
		IHybridDescriptor lipo = new LipoTopologicDescriptor();
		value = refract.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = refract.getAtomContainer();
		value = elect.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = elect.getAtomContainer();
		value = lipo.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = lipo.getAtomContainer();
		this.graph = aux;
	}

	public void calcularIndicesTopograficos() throws CloneNotSupportedException, CDKException {
		DescriptorValue value;
		IAtomContainer aux = (IAtomContainer) this.graph.clone();
		IHybridDescriptor refract = new RefractedTopographicDescriptor();
		IHybridDescriptor elect = new ElectroTopographicDescriptor();
		IHybridDescriptor lipo = new LipoTopographicDescriptor();
		value = refract.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = refract.getAtomContainer();
		value = elect.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = elect.getAtomContainer();
		value = lipo.calculate(aux);
		if (value.getException() != null)
			throw new CDKException(value.getException().getMessage());
		aux = lipo.getAtomContainer();
		this.graph = aux;
	}

	public double sumaTotalDescriptorTopographical(String index) {
		double suma = 0;
		for (int i = 0; i < graph.getAtomCount(); i++) {
			if (!graph.getAtom(i).getSymbol().equals("H")) {
				if(index == "E")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
				else if(index == "R")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
				else if(index == "L")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
				else if(index == "ER") 
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC), 2));
				else if(index == "EL") 
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC), 2));
				else if(index == "RL")
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC), 2));		
				else
				   suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC), 2) +
						   Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC), 2) + 
						   Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC), 2));
			}
		}
		return suma;
	}
	
	public double sumaTotalDescriptorTopographicalNorm(String index) {
		double suma = 0;
		for (int i = 0; i < graph.getAtomCount(); i++) {
			if (!graph.getAtom(i).getSymbol().equals("H")) {
				if(index == "E")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM);
				else if(index == "R")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM);
				else if(index == "L")
				   suma += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM);
				else if(index == "ER") 
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2));
				else if(index == "EL") 
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
				else if(index == "RL")
					suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) +
							Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));		
				else
				   suma += Math.sqrt(Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM), 2) +
						   Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) + 
						   Math.pow((Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
			}
		}
		return suma;
	}
	
	public double[] sumaTotalDescriptorTopological() {
		double[] suma = { 0, 0, 0 };
		for (int i = 0; i < graph.getAtomCount(); i++) {
			if (!graph.getAtom(i).getSymbol().equals("H")) {
				suma[0] += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOLOGICAL);
				suma[1] += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOLOGICAL);
				suma[2] += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOLOGICAL);
			}
		}
		return suma;
	}

	private double valoresRefractoCarbonos() {
		double valor = 0;
		for (int i = 0; i < graph.getAtomCount(); i++) 
			if (graph.getAtom(i).getSymbol().equals("C"))
				valor += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
		return valor;
	}

	private double valorElectroHeteroatomos() {
		double valor = 0;
		for (int i = 0; i < graph.getAtomCount(); i++)
			if (!graph.getAtom(i).getSymbol().equals("H") && !graph.getAtom(i).getSymbol().equals("C"))
				valor += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
		return valor;
	}

	private double valorLipoHeteroatomos() {
		double valor = 0;
		for (int i = 0; i < graph.getAtomCount(); i++)
			if (!graph.getAtom(i).getSymbol().equals("H") && !graph.getAtom(i).getSymbol().equals("C"))
				valor += (Double) graph.getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
		return valor;
	}

	private double valoresHidrogenosHeteroatomos() {
		double valor = 0;
		for (int i = 0; i < graph.getAtomCount(); i++)
			if (!graph.getAtom(i).getSymbol().equals("H") && !graph.getAtom(i).getSymbol().equals("C"))
				valor += graph.getAtom(i).getImplicitHydrogenCount() * 0.23;
		return valor;
	}

	public double getPesoETPG() {
		return valorElectroHeteroatomos() / sumaTotalDescriptorTopographical("E");
	}

	public double getPesoRTPG() {
		return valoresRefractoCarbonos() / sumaTotalDescriptorTopographical("R");
	}

	public double getPesoLTPG() {
		return (valorLipoHeteroatomos() + valoresHidrogenosHeteroatomos()) / sumaTotalDescriptorTopographical("L");
	}
	
	public double getPesoETLG() {
		return valorElectroHeteroatomos() / sumaTotalDescriptorTopological()[0];
	}

	public double getPesoRTLG() {
		return valoresRefractoCarbonos() / sumaTotalDescriptorTopological()[1];
	}

	public double getPesoLTLG() {
		return (valorLipoHeteroatomos() + valoresHidrogenosHeteroatomos()) / sumaTotalDescriptorTopological()[2];
	}

	public int getCantAtmosPesados() {
		int cant = 0;
		for (int i = 0; i < graph.getAtomCount(); i++)
			if (!graph.getAtom(i).getSymbol().equals("H"))
				cant++;
		return cant;
	}
	
	public ArrayList<IAtom> getAtmosPesados() {
		ArrayList<IAtom> atoms = new ArrayList<IAtom>();
		for (int i = 0; i < graph.getAtomCount(); i++)
			if (!graph.getAtom(i).getSymbol().equals("H"))
				atoms.add(graph.getAtom(i));
		return atoms;
	}

	public int getCantAtmos() {
		return graph.getAtomCount();
	}

	public IAtomContainer getGraph() {
		return graph;
	}

	public int getCantDescriptorCenter() {
		return reducedGraph.getNumeroCentrosDescriptores();
	}

	@SuppressWarnings("rawtypes")
	public Set getListDescriptorCenter() {
		return reducedGraph.getAllCentrosDescriptores();
	}

	public double getDistancia(int pos1, int pos2) {
		return distancia[pos1][pos2];
	}

	// *****************************************************************************************************
	// ************************IMPLEMENTADO POR JUAN LUIS
	// PANEQUE*******************************************
	// *******************************************************************************************************

	/**
	 * Devuelve los fragmentos moleculardes de tamaño M centros de descriptores.
	 * 
	 * @param M
	 * @return
	 */
	public ArrayList<FragmentoMolecular> getFragmentosMolecularesOrdenN(int M) {
		ArrayList<FragmentoMolecular> fragmentos = new ArrayList<FragmentoMolecular>();
		Object[] listDescriptorCenter = reducedGraph.getAllCentrosDescriptores().toArray();
		int N = listDescriptorCenter.length;
		for (int n = 0; n < 1 << N; n++) {
			int v = n - ((n >> 1) & 0x55555555); // Cuenta el número de bits a 1 // (sólo para 32 bits)							
			v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
			v = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
			FragmentoMolecular fragmento = new FragmentoMolecular();
			ArrayList<Integer> centros = new ArrayList<Integer>();
			if (v == M) {
				for (int i = 0; i < N; i++)
					if ((n & (1 << i)) != 0) {
						fragmento.addCentroDescriptor((CentroDescriptor) listDescriptorCenter[i]);
						centros.add(i);
					}
				fragmento.setMolecula(this);
				fragmentos.add(fragmento);
			}
		}
		
		for (int i=0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				FragmentoMolecular fragmento = new FragmentoMolecular();
				boolean tt = includeCentroDescriptor((CentroDescriptor) listDescriptorCenter[i], (CentroDescriptor) listDescriptorCenter[j]);
				if (!tt) {
					fragmento.addCentroDescriptor((CentroDescriptor) listDescriptorCenter[i]);
					fragmento.addCentroDescriptor((CentroDescriptor) listDescriptorCenter[j]);
					fragmento.setMolecula(this);
				    fragmentos.add(fragmento);	
				}
				
			}
		    
		}
			
				
		
		return fragmentos;
	}

	/**
	 * Devuelve el fragmento molecular correspondiente a una secuencia de atomos
	 * 
	 * @param secuence
	 * @return
	 */
	public FragmentoMolecular getFragmentoOfSec(String secuence) {
		ArrayList<FragmentoMolecular> fragmentos = getFragmentosMolecularesOrdenN(2);
		for (FragmentoMolecular fragmento : fragmentos) 
			if (fragmento.isSecuenceCorrect(secuence))
				return fragmento;
		return null;
	}

	/**
	 * Devuelve los fragmentos moleculares similares al fragmento base introducido utilizando una funcion especificada
	 * 
	 * @param type
	 * @param max
	 * @param error
	 * @param fragmentoBase
	 * @return similares
	 * @throws CDKException 
	 */
	public ArrayList<FragmentoMolecular> getSimilarFragments(TypeSimilaryFunction type, boolean max, double error, FragmentoMolecular fragmentoBase, char tipo, char indice) throws CDKException {
		ArrayList<FragmentoMolecular> fragmentos = getFragmentosMolecularesOrdenN(2);
		ArrayList<FragmentoMolecular> similares = new ArrayList<FragmentoMolecular>();
		double index1 = 0.0D;
		double index2 = 0.0D;
		double index = 0.0D;
		double[] frag1 = null;
		double[] frag2 = null;
		for (FragmentoMolecular fragment : fragmentos) {
			if(tipo=='P'){
			  frag1 = fragmentoBase.vectorValoresTopograficos(indice);
			  frag2 = fragment.vectorValoresTopograficos(indice); 
		    } else {
		      frag1 = fragmentoBase.vectorValoresTopologicos(indice);
			  frag2 = fragment.vectorValoresTopologicos(indice);
		    }
			
			double[] ffrag1 = {frag1[0], frag1[1]+frag1[4], frag1[2]+frag1[5], frag1[3]+frag1[6]};
			double[] ffrag2 = {frag2[0], frag2[1]+frag2[4], frag2[2]+frag2[5], frag2[3]+frag2[6]};
			
			index1 = similaryFunction.CalculeSimilaryFunction(type, ffrag1, ffrag2);
			/*double[] frag3;
			if (frag2.length == 7)
				frag3 = new double[] { frag2[0], frag2[3], frag2[4], frag2[5], frag2[1], frag2[2], frag2[6] };
			else
				frag3 = new double[] { frag2[0], frag2[2], frag2[1] };
			index2 = similaryFunction.CalculeSimilaryFunction(type, frag1, frag3);
			if (max) {
				index = Math.max(index1, index2);
				if (index > error) {
					fragment.setIndexSimilary(index);
					similares.add(fragment);
				}
			} else {
				index = Math.min(index1, index2);
				if (index < error) {
					fragment.setIndexSimilary(index);
					similares.add(fragment);
				}
			}*/
			
			if (index1 > error) {
				fragment.setIndexSimilary(index1);
				similares.add(fragment);
			}
			
		}
		return similares;
	}
	
	/*public AtomsPMC getAtomsPMC(TypeSimilaryFunction funcSim, boolean max, double error, Molecule mol, char indice, char metodo) throws CDKException {
		ArrayList<IAtom> listAtoms = new ArrayList<IAtom>();
		for(IAtom atom : getGraph().atoms())
			if(!atom.getSymbol().equals("H"))
			  listAtoms.add(atom);
		
		ArrayList<IAtom> listAtomsComp = mol.getAtomsSelected();
		for(IAtom atom : mol.getGraph().atoms())
			if(!atom.getSymbol().equals("H"))
			  listAtomsComp.add(atom);
		
		int filas = listAtoms.size();
		int columns = listAtomsComp.size();
		MatrixCompared matrix = new MatrixCompared(filas, columns);
		for (int i = 0; i < filas; i++) {
			double [] listProperty1 = new double [3];
			listProperty1[0] = (Double) listAtoms.get(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
			listProperty1[1] = (Double) listAtoms.get(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
			listProperty1[2] = (Double) listAtoms.get(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
			for (int j = 0; j < columns; j++) {
				double [] listProperty2 = new double [3];
				listProperty2[0] = (Double) listAtomsComp.get(j).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC);
				listProperty2[1] = (Double) listAtomsComp.get(j).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC);
				listProperty2[2] = (Double) listAtomsComp.get(j).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC);
				double index = similaryFunction.CalculeSimilaryFunction(funcSim, listProperty1, listProperty2);
				matrix.addData(new Data(index, i, j));
			}
		}
		
		String directorio = "D:\\Doctorado";
		
		ExcelStadistics excel = new ExcelStadistics();
		try {
			excel.createMatrixSimilariry(listAtoms, listAtomsComp, matrix, directorio);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ArrayList<Data> listCompared = matrix.Compared(max, error);
		
		ArrayList<IAtom> listMol1 = new ArrayList<IAtom>();
		ArrayList<IAtom> listMol2 = new ArrayList<IAtom>();
		
		/*FragmentoMolecular fragm1 = new FragmentoMolecular();
		fragm1.setMolecula(this);
		FragmentoMolecular fragm2 = new FragmentoMolecular();
		fragm2.setMolecula(mol);
		for (Data d : listCompared) {
			listMol1.add(listAtoms.get(d.getPosX(0)));
			listMol2.add(listAtomsComp.get(d.getPosY(0)));
		}
		
		/*int filas1 = fragm1.getGradoFragmento();
		int columns1 = fragm1.getGradoFragmento();
		MatrixCompared mol1 = new MatrixCompared(filas1, columns1);
		MatrixCompared mol2 = new MatrixCompared(filas1, columns1);
		MatrixCompared molDist = new MatrixCompared(filas1, columns1);
		
		
		ArrayList<CentroDescriptor> listFragm1 = fragm1.getCentDescriptores();
		ArrayList<CentroDescriptor> listFragm2 = fragm2.getCentDescriptores();
		
		FragmentoMolecular Fragm1 = new FragmentoMolecular();
		FragmentoMolecular Fragm2 = new FragmentoMolecular();
		
		for (CentroDescriptor e : listFragm1) {
			Fragm1.addCentroDescriptor(e);
		}
		
		for (CentroDescriptor e : listFragm2) {
			Fragm2.addCentroDescriptor(e);
		}
		
		AtomsPMC fpmc = new AtomsPMC (listMol1, listMol2);
		
		
		
		for (int i = 0; i < listFragm1.size(); i++) {
			for (int j = 0; j < listFragm1.size(); j++) {
				double index1 = listFragm1.get(i).getCentroMasa().distance(listFragm1.get(j).getCentroMasa());
				mol1.addData(new Data(index1, i, j, false));
				double index2 = listFragm2.get(i).getCentroMasa().distance(listFragm2.get(j).getCentroMasa());
				mol2.addData(new Data(index2, i, j, false));
			}
		}
			
		molDist = mol1.Similary(funcSim, mol2);
		
		ArrayList<Data> listDistCompared = molDist.ComparedMatix(max, error);
		
		if (listDistCompared.size() > 0) {
			FragmentoMolecular listFragm1Final = new FragmentoMolecular();
			FragmentoMolecular listFragm2Final = new FragmentoMolecular();
			for (Data d : listDistCompared) {
				listFragm1Final.addCentroDescriptor(listFragm1.get(d.getPosX()));
				listFragm1Final.addCentroDescriptor(listFragm1.get(d.getPosY()));
				listFragm2Final.addCentroDescriptor(listFragm2.get(d.getPosX()));
				listFragm2Final.addCentroDescriptor(listFragm2.get(d.getPosY()));
			}
			
			FragmentosPMC fpmc = new FragmentosPMC (listFragm1Final, listFragm2Final);
			
			//fpmc.setIndexSimilary(index);
			return fpmc;
		}
		
		return fpmc;
	}*/
	
	public FragmentosPMC searchFragmentPMCCD(ArrayList<CentroDescriptor> listCentroDescriptors, TypeSimilaryFunction funcSim, double error, Molecule mol, String indice) throws CDKException{
		FragmentosPMC fpmc = null;
		Object[] centroDescriptores = listCentroDescriptors.toArray();
		Object[] centroDescriptoresComp = mol.getListDescriptorCenter().toArray();
		int filas = centroDescriptores.length;
		int columns = centroDescriptoresComp.length;
		MatrixCompared matrix = new MatrixCompared(filas, columns);
		for (int i = 0; i < filas; i++)
			for (int j = 0; j < columns; j++) {
				double index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, ((CentroDescriptor) centroDescriptores[i]).getValorTotalVector(indice), ((CentroDescriptor) centroDescriptoresComp[j]).getValorTotalVector(indice));
				ArrayList<Double> value = new ArrayList<Double>();
				value.add(index);
				ArrayList<Integer> x = new ArrayList<Integer>();
				x.add(i);
				ArrayList<Integer> y = new ArrayList<Integer>();
				y.add(j);
				matrix.addData(new Data(value, x, y));
			}
        ArrayList<Data> listCompared = matrix.ComparedCD(error, centroDescriptores, centroDescriptoresComp, indice);
		
		if(!listCompared.isEmpty()){
			FragmentoMolecular fragm1 = new FragmentoMolecular();
			fragm1.setMolecula(this);
			FragmentoMolecular fragm2 = new FragmentoMolecular();
			fragm2.setMolecula(mol);
			String cadenaFrag1 = "Frag1: ";
			String cadenaFrag2 = "Frag2: ";
			for (Data d : listCompared) {
				if(!fragm1.ContainCenterDescriptor((CentroDescriptor) centroDescriptores[d.getPosX(0)]))
				   fragm1.addCentroDescriptor((CentroDescriptor) centroDescriptores[d.getPosX(0)]);
		
				if(!fragm2.ContainCenterDescriptor((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]))
				   fragm2.addCentroDescriptor((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]);
				
			    cadenaFrag1 += ((CentroDescriptor) centroDescriptores[d.getPosX(0)]).getNombreCD() + "-";
			    cadenaFrag2 += ((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]).getNombreCD() + "-";
			}
			
			System.out.println(cadenaFrag1);
			System.out.println(cadenaFrag2);
			
			if(fragm1.getGradoFragmento() !=0 || fragm2.getGradoFragmento() !=0){
			  fpmc = new FragmentosPMC (fragm1, fragm2);
			}
		}
		return fpmc;
	}
	
	public FragmentosPMC getFragmentoPMCCD(TypeSimilaryFunction funcSim, double error, Molecule mol, String indice) throws CDKException {
		FragmentosPMC fpmc = null;
		Object[] centroDescriptores = getListDescriptorCenter().toArray();
		Object[] centroDescriptoresComp = mol.getListDescriptorCenter().toArray();
		int filas = centroDescriptores.length;
		int columns = centroDescriptoresComp.length;
		MatrixCompared matrix = new MatrixCompared(filas, columns);
		for (int i = 0; i < filas; i++) {
			for (int j = 0; j < columns; j++) {
				//double indexE = similaryFunction.CalculeSimilaryFunction(funcSim, ((CentroDescriptor) centroDescriptores[i]).getValoresTopogTotal('E'), ((CentroDescriptor) centroDescriptoresComp[j]).getValoresTopogTotal('E'));
				//double indexR = similaryFunction.CalculeSimilaryFunction(funcSim, ((CentroDescriptor) centroDescriptores[i]).getValoresTopogTotal('R'), ((CentroDescriptor) centroDescriptoresComp[j]).getValoresTopogTotal('R'));
				//double indexL = similaryFunction.CalculeSimilaryFunction(funcSim, ((CentroDescriptor) centroDescriptores[i]).getValoresTopogTotal('L'), ((CentroDescriptor) centroDescriptoresComp[j]).getValoresTopogTotal('L'));
				double index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, ((CentroDescriptor) centroDescriptores[i]).getValorTotalVector(indice), ((CentroDescriptor) centroDescriptoresComp[j]).getValorTotalVector(indice));
				//System.out.println(((CentroDescriptor) centroDescriptores[i]).getNombreCD() + " - " + ((CentroDescriptor) centroDescriptores[i]).getValorTotalVector() + " - " + ((CentroDescriptor) centroDescriptoresComp[j]).getNombreCD() + " - " + ((CentroDescriptor) centroDescriptoresComp[j]).getValorTotalVector() + " - " + index);	
				ArrayList<Double> value = new ArrayList<Double>();
				value.add(index);
				ArrayList<Integer> x = new ArrayList<Integer>();
				x.add(i);
				ArrayList<Integer> y = new ArrayList<Integer>();
				y.add(j);
				matrix.addData(new Data(value, x, y));
			}
		}
		
	/*	String directorio = "D:\\Doctorado";
		
		ExcelStadistics excel = new ExcelStadistics();
		try {
			excel.createMatrixSimilariry(centroDescriptores, centroDescriptoresComp, matrix, directorio);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
		
		ArrayList<Data> listCompared = matrix.ComparedCD(error, centroDescriptores, centroDescriptoresComp, indice);
		
		if(!listCompared.isEmpty()){
			FragmentoMolecular fragm1 = new FragmentoMolecular();
			fragm1.setMolecula(this);
			FragmentoMolecular fragm2 = new FragmentoMolecular();
			fragm2.setMolecula(mol);
			String cadenaFrag1 = "Frag1: ";
			String cadenaFrag2 = "Frag2: ";
			for (Data d : listCompared) {
				if(!fragm1.ContainCenterDescriptor((CentroDescriptor) centroDescriptores[d.getPosX(0)]))
				   fragm1.addCentroDescriptor((CentroDescriptor) centroDescriptores[d.getPosX(0)]);
		
				if(!fragm2.ContainCenterDescriptor((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]))
				   fragm2.addCentroDescriptor((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]);
				
			    cadenaFrag1 += ((CentroDescriptor) centroDescriptores[d.getPosX(0)]).getNombreCD() + "-";
			    cadenaFrag2 += ((CentroDescriptor) centroDescriptoresComp[d.getPosY(0)]).getNombreCD() + "-";
			}
			
			System.out.println(cadenaFrag1);
			System.out.println(cadenaFrag2);
			
			if(fragm1.getGradoFragmento() !=0 || fragm2.getGradoFragmento() !=0){
			  fpmc = new FragmentosPMC (fragm1, fragm2);
			}
		}
		
		return fpmc;
	}
	
	public FragmentosPMC getFragmentoPMCAtoms(TypeSimilaryFunction funcSim, boolean max, double error, Molecule mol, String indice) throws CDKException {
		FragmentosPMC fpmc = null;
		Object[] atoms = getAtmosPesados().toArray();
		Object[] atomsComp = mol.getAtmosPesados().toArray();
		int filas = atoms.length;
		int columns = atomsComp.length;
		MatrixCompared matrix = new MatrixCompared(filas, columns);
		for (int i = 0; i < filas; i++) {
			for (int j = 0; j < columns; j++) {
				double index =0.0;
				if(indice.equals("E"))
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, (Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM), (Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM));
				else if (indice.equals("R"))
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, (Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), (Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM));
				else if (indice.equals("L"))
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, (Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), (Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM));
				else if (indice.equals("EL")){
					double valor1 = Math.sqrt(Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2)+ Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
					double valor2 = Math.sqrt(Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2)+ Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, valor1, valor2);
				}else if (indice.equals("ER")){
					double valor1 = Math.sqrt(Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2));
					double valor2 = Math.sqrt(Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2));
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, valor1, valor2);
				}else if (indice.equals("RL")){
					double valor1 = Math.sqrt(Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
					double valor2 = Math.sqrt(Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM), 2));
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, valor1, valor2);
				}else{
					double valor1 = Math.sqrt(Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) + 
							        Math.pow((Double)((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM),2));
					double valor2 = Math.sqrt(Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM),2) + Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM), 2) +
							        Math.pow((Double)((IAtom) atomsComp[j]).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM),2));
					index = similaryFunction.CalculeSimilaryFunction(TypeSimilaryFunction.TANIMOTO, valor1, valor2);
				}
				
				/*System.out.println(((IAtom) atoms[i]).getSymbol()+ ((IAtom) atoms[i]).getID() + " - " + ((IAtom) atoms[i]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + 
						           " - " + ((IAtom) atomsComp[j]).getSymbol()+ ((IAtom) atomsComp[j]).getID() + " - " + ((IAtom) atoms[j]).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + 
						           " - " + index);*/	
				ArrayList<Double> value = new ArrayList<Double>();
				value.add(index);
				ArrayList<Integer> x = new ArrayList<Integer>();
				x.add(i);
				ArrayList<Integer> y = new ArrayList<Integer>();
				y.add(j);
				matrix.addData(new Data(value, x, y));
			}
		}
		
	   String directorio = "D:\\Doctorado";
		
		ExcelStadistics excel = new ExcelStadistics();
		try {
			excel.createMatrixSimilariryAtoms(atoms, atomsComp, matrix, directorio);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ArrayList<Data> listCompared = matrix.ComparedAtoms(max, error, atoms, atomsComp, indice);
		
		if(!listCompared.isEmpty()){
			FragmentoMolecular fragm1 = new FragmentoMolecular();
			fragm1.setMolecula(this);
			FragmentoMolecular fragm2 = new FragmentoMolecular();
			fragm2.setMolecula(mol);
			String cadenaFrag1 = "Frag1: ";
			String cadenaFrag2 = "Frag2: ";
			for (Data d : listCompared) {
				if(!fragm1.ContainAtom((IAtom)atoms[d.getPosX(0)]))
				   fragm1.addAtoms((IAtom)atoms[d.getPosX(0)]);
		
				if(!fragm2.ContainAtom((IAtom)atomsComp[d.getPosY(0)]))
				   fragm2.addAtoms((IAtom)atomsComp[d.getPosY(0)]);
				
			    cadenaFrag1 += ((IAtom) atoms[d.getPosX(0)]).getSymbol() + ((IAtom) atoms[d.getPosX(0)]).getID() + "-";
			    cadenaFrag2 += ((IAtom) atomsComp[d.getPosY(0)]).getSymbol() + ((IAtom) atomsComp[d.getPosY(0)]).getID() + "-";
			}
			
			System.out.println(cadenaFrag1);
			System.out.println(cadenaFrag2);
			
			if(fragm1.getCantAtoms() !=0 || fragm2.getCantAtoms() !=0){
			  fpmc = new FragmentosPMC (fragm1, fragm2);
			}
		}
		
		return fpmc;
	}
	
	
	public void imprmirMatrix(MatrixCompared m){
		for(int i=0; i< m.getFile(); i++){
			String cadena = "";
			for(int j=0; j< m.getColumn(); j++)
				cadena += m.getData()[i][j].getValue() + ";";
			//System.out.println(cadena);
		}
			
	}
	
	public void clearListAtomSelected(){
		listAtomSelected.clear();;
	}
	
	public boolean isAtomSelected(IAtom atom){
		boolean flag = false;
		for(IAtom a : listAtomSelected)
			if(a.getSymbol().equals(atom.getSymbol()) && a.getID().equals(atom.getID())){
				flag = true;
				break;
			}	
		return flag;
	}
	
	

	public ArrayList<Data> getListCDSimilares(Object[] centDescComp, char indice) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public ArrayList<IAtom> getAtoms(){
		ArrayList<IAtom> listAtoms =  new ArrayList<IAtom>();
		for(IAtom atom : getGraph().atoms())
			listAtoms.add(atom);
		return listAtoms;
	}
	
	public void normalizationDescriptor(Molecule molecule){
		
		ArrayList<Double> list = new ArrayList<Double>();
		for (int i=0; i< molecule.getCantAtmos(); i++) {
			if (!molecule.getGraph().getAtom(i).getSymbol().equals("H")) {
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC));
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC));
				list.add((Double) molecule.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC));
			}
		}
		
		double mayorMol2 = -0.00009;
		for (int i = 0; i < list.size(); i++)
			if (list.get(i)<0 && list.get(i) < mayorMol2)
				mayorMol2 = list.get(i);
		
		ArrayList<Double> list1 = new ArrayList<Double>();
		for (int i=0; i< this.getCantAtmos(); i++) {
			if (!this.getGraph().getAtom(i).getSymbol().equals("H")) {
				list1.add((Double) this.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC));
				list1.add((Double) this.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC));
				list1.add((Double) this.getGraph().getAtom(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC));
			}
		}
		
		double mayorMol1 = -0.00009;
		for (int i = 0; i < list1.size(); i++)
			if (list1.get(i)<0 && list1.get(i) < mayorMol1)
				mayorMol1 = list1.get(i);
			    
		if(mayorMol1 < mayorMol2){
			for(int i=0; i<this.getCantAtmos(); i++)
				if(!this.getAtoms().get(i).getSymbol().equals("H")){
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) + Math.abs(mayorMol1) + 1);
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + Math.abs(mayorMol1) + 1);
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) + Math.abs(mayorMol1) + 1);
				}
			for(int i=0; i<molecule.getCantAtmos(); i++)
				if(!molecule.getAtoms().get(i).getSymbol().equals("H")){
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) + Math.abs(mayorMol1) +1);
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + Math.abs(mayorMol1) +1);
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) + Math.abs(mayorMol1) +1);
				}
		} else {
			for(int i=0; i<this.getCantAtmos(); i++)
				if(!this.getAtoms().get(i).getSymbol().equals("H")){
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
					this.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM, (Double) this.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
				}
			for(int i=0; i<molecule.getCantAtmos(); i++)
				if(!molecule.getAtoms().get(i).getSymbol().equals("H")){
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.LIPOTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
					molecule.getAtoms().get(i).setProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHICNORM, (Double) molecule.getAtoms().get(i).getProperty(IHybridDescriptor.Descriptor.REFRACTEDTOPOGRAPHIC) + Math.abs(mayorMol2) +1);
				}
		}
	
	}
}
