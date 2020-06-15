package domain;

import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import org._3pq.jgrapht.edge.UndirectedWeightedEdge;
import org._3pq.jgrapht.graph.SimpleGraph;

public class ReducedGraph {

	private SimpleGraph moleculaCD;
	private int nAnillos,nCluster4,nCluster3,nMetilo,nMetileno,nMetino,nHet, nCConec, id;

	public ReducedGraph() {
		moleculaCD = new SimpleGraph();
		nAnillos = nCluster4 = nCluster3 = nMetilo = nMetileno = nMetino= nHet = nCConec = id = 0;		
	}

	public void addCentroDescriptor(CentroDescriptor cd){
		cd.setId(id++);
		moleculaCD.addVertex(cd);
		actualizarCantidades(cd.getTipo());
	}

	public void addVariosCentrosDescriptores(List<CentroDescriptor> c){
		@SuppressWarnings("rawtypes")
		Iterator a = c.iterator();
		while (a.hasNext()) {
			CentroDescriptor cd = (CentroDescriptor)a.next();
			addCentroDescriptor(cd);
		}
	}

	private void actualizarCantidades(TipoCentroDescriptor tipo) {
		switch (tipo) {
		case ANILLO: nAnillos++; break;
		case CLUSTER4: nCluster4++; break;
		case CLUSTER3: nCluster3++; break;
		case METILENO: nMetileno++; break;
		case METILO: nMetilo++; break;
		case METINO: nMetino++; break;
		case HETEROATOMO: nHet++; break;
		case CARBONOCONECTADO: nCConec++; break;
		default:
			break;		
		}		
	}

	@SuppressWarnings("rawtypes")
	public Set getAllCentrosDescriptores(){
		return moleculaCD.vertexSet();
	}	
	

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public Set getCentDescripPorTipo(TipoCentroDescriptor tipo){
		Iterator cd = moleculaCD.vertexSet().iterator();
		Set nodos = new LinkedHashSet();
		while (cd.hasNext()) {
			CentroDescriptor e = (CentroDescriptor) cd.next();
			if(e.getTipo().equals(tipo))
				nodos.add(e);
		}

		return nodos;
	}
	
	public int getNumeroCentrosDescriptores(){
		return moleculaCD.vertexSet().size();
	}

	/*
	 * Grafo reducido topogr�fico
	 */
	public void addCaminoTopografico(CentroDescriptor cd1, CentroDescriptor cd2, double dist){
		UndirectedWeightedEdge edge = new UndirectedWeightedEdge(cd1, cd2, dist);
		moleculaCD.addEdge(edge);
	}

	public double getDistancia(CentroDescriptor cd1, CentroDescriptor cd2){
		return moleculaCD.getEdge(cd1, cd2).getWeight();
	}
	
	public String caracteristicasGenerales(){
		StringBuilder sb = new StringBuilder();
		/*sb.append("Anillos " + nAnillos + "\n");
		sb.append("Cluster4 " + nCluster4 + "\n");
		sb.append("Cluster3 " + nCluster3 + "\n");
		sb.append("Heteroatomos " + nHet + "\n");
		sb.append("Metino " + nMetino + "\n");
		sb.append("Metileno " + nMetileno + "\n");
		sb.append("Metilo " + nMetilo + "\n");*/
		int total = nAnillos + nCluster4 + nCluster3 + nHet + nMetino + nMetileno + nMetilo + nCConec;
		sb.append("Total " + total + "\n");
		sb.trimToSize();
		return sb.toString();
	}
	
	/*private StringBuilder caracteristicasCentrosDescp(String molecula){
		StringBuilder stb = new StringBuilder();
		Iterator tempIterator = new DepthFirstIterator(moleculaCD);
	    while(tempIterator.hasNext()){
	     CentroDescriptor aux = (CentroDescriptor)tempIterator.next();
	     String s = molecula + ";";
	    	 s += aux.toStringLeft();
	     stb.append(s);
	     stb.append("\n");
	    }
		stb.trimToSize();
		return stb;
	}*/
	
/*	private StringBuilder caracteristicasEnlaces(String molecula){
		StringBuilder stb = new StringBuilder();
		int n = moleculaCD.vertexSet().size();
		Object[] cd = moleculaCD.vertexSet().toArray();
		double suma[]= {0, 0, 0};
		for (int i = 0; i < n; i++) {
			CentroDescriptor temp = (CentroDescriptor) cd[i];
			for (int j = 0; j < temp.getCantidadAtomos(); j++) {
				System.out.println(temp.getSimbolAtomo(j) + " " +temp.getIdAtomo(j));
				if(!temp.getSimbolAtomo(j).equals("H"))
				{
					stb.append(molecula + ";");
					stb.append(temp.getSimbolAtomo(j) + ";");
					stb.append(temp.getIdAtomo(j) + ";");
					stb.append(temp.getValoresTopogAtomo(j)[0] + ";");
					stb.append(temp.getValoresTopogAtomo(j)[1] + ";");
					stb.append(temp.getValoresTopogAtomo(j)[2] + "\n");
				}
				
			}
		}
		
		for(int i = 0; i < n - 1; i++){
			for(int j = i+1; j < n; j++){
				CentroDescriptor cd1 = (CentroDescriptor) cd[i];
				CentroDescriptor cd2 = (CentroDescriptor) cd[j];
				double dist1 = moleculaCD.getEdge(cd1, cd2).getWeight();
				stb.append(molecula + ";");
				stb.append(cd1.getNombreCD() + "-" + cd2.getNombreCD() + ";");
				stb.append(cd1.getSecuencia() + "-" + cd2.getSecuencia() + ";");
				stb.append(cd1.toStringRight() + ";");
				stb.append(dist1 + ";");
				stb.append(cd2.toStringRight() + ";");
				stb.append("\n");
			}
		}
		for (int j = 0; j < cd.length ; j++) {
			suma[0] += ((CentroDescriptor) cd[j]).valoresTopogTotal[0];
			suma[1] += ((CentroDescriptor) cd[j]).valoresTopogTotal[1];
			suma[2] += ((CentroDescriptor) cd[j]).valoresTopogTotal[2];
		}
		stb.append(molecula + ";");
		stb.append(suma[0] + ";");
		stb.append(suma[1] + ";");
		stb.append(suma[2] + ";");
		stb.append("\n");
	    stb.trimToSize();
		return stb;
	}*/

	/*public String toString(String molecula) {
		StringBuilder builder = new StringBuilder();
		//builder.append(caracteristicasGenerales());
		//builder.append(caracteristicasCentrosDescp(molecula));
		builder.append(caracteristicasEnlaces(molecula.substring(0, molecula.length()-4)));
		builder.trimToSize();
		return builder.toString();
	}*/
	
	
	

}
