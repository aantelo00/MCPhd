package util;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

import org.openscience.cdk.interfaces.IAtom;

import domain.CentroDescriptor;

public class MatrixCompared {

	private Data[][] matrix;
	private int file;
	private int column;
	private double similitudPromedio;

	public MatrixCompared(int file, int column) {
		matrix = new Data[file][column];
		this.file = file;
		this.column = column;
	}

	public void addData(Data data) {
		matrix[data.getPosX(0)][data.getPosY(0)] = data;
	}
	
	public Data[][] getData(){
		return matrix;
	}
	
	public Data getData(int x, int y){
      return matrix[x][y];	
   }
	
	public int getFile(){
		return file;
	}
	
	public int getColumn(){
		return column;
	}
	
	public Data getMaximo(){
		double mayor = Double.MIN_VALUE;
		int posX = -1, posY = -1;
		for (int i = 0; i < file; i++) 
			for (int j = 0; j < column; j++) {
				if (matrix[i][j].getValue().get(0) >= mayor && matrix[i][j].getFlag()==false){
					mayor = matrix[i][j].getValue().get(0);
					posX = i;
					posY = j;
				}
			}
		matrix[posX][posY].setFlag(true);
		return matrix[posX][posY];
	}

	/**
	 * @param value
	 * @return La lista de centros de descriptores cuya semejanza sea mayor/menor que el valor especificado
	 * 
	 */
	
	public static double round(double value, int places){
		if(places < 0) throw new IllegalArgumentException();
		
		BigDecimal bd = new BigDecimal(value);
		bd = bd.setScale(places, RoundingMode.HALF_UP);
		return bd.doubleValue();
	}
	
	public ArrayList<Data> ComparedCD(double error, Object[] centroDescriptores, Object[] centroDescriptoresComp, String indice) {
		ArrayList<Data> lista = new ArrayList<Data>();
		boolean flag = true;
		while(flag){
			Data dataTemp = getMaximo();
			if(round(dataTemp.getValue().get(0),2)>=error )
				lista.add(dataTemp);
			else
				flag = false;
		}
		
        ArrayList<ArrayList<Data>> listaFinal = new ArrayList<ArrayList<Data>>();
		//String cadena1 ="";
        for(int i=0; i<lista.size(); i++){
        	ArrayList<Data> listaTemp = new ArrayList<Data>();
			listaTemp.add(lista.get(i));
			 //cadena1 = ((CentroDescriptor) centroDescriptores[lista.get(i).getPosX(0)]).getNombreCD()+ "("+((CentroDescriptor) centroDescriptores[lista.get(i).getPosX(0)]).getSecuencia() +")" + "-" +
			 //		  ((CentroDescriptor) centroDescriptoresComp[lista.get(i).getPosY(0)]).getNombreCD() + "(" + ((CentroDescriptor) centroDescriptoresComp[lista.get(i).getPosY(0)]).getSecuencia() + ")"+ " ";
			double index =0;
			for(int j=0; j<lista.size(); j++){
				if(i!=j){
				   boolean flagTemp = true;
				   for(int k=0; k<listaTemp.size(); k++){
					  double valor1 = ((CentroDescriptor) centroDescriptores[lista.get(j).getPosX(0)]).getCentroMasa().distance(((CentroDescriptor) centroDescriptores[listaTemp.get(k).getPosX(0)]).getCentroMasa());
					  double valor2 = ((CentroDescriptor) centroDescriptoresComp[lista.get(j).getPosY(0)]).getCentroMasa().distance(((CentroDescriptor) centroDescriptoresComp[listaTemp.get(k).getPosY(0)]).getCentroMasa());
					
					  index = 1-Math.abs(valor1-valor2)/(valor1+valor2);
					  if(index < 0.85)
						  flagTemp = false; 
				   }
				   if(flagTemp){
					   listaTemp.add(lista.get(j));
					   //cadena1 += ((CentroDescriptor) centroDescriptores[lista.get(j).getPosX(0)]).getNombreCD()+ "("+((CentroDescriptor) centroDescriptores[lista.get(j).getPosX(0)]).getSecuencia() +")"+"-" +
					   //			  ((CentroDescriptor) centroDescriptoresComp[lista.get(j).getPosY(0)]).getNombreCD()+ "(" + ((CentroDescriptor) centroDescriptoresComp[lista.get(j).getPosY(0)]).getSecuencia(); 
				   }	   
			    }
		   }
			listaFinal.add(listaTemp);
		}
		
		if(!listaFinal.isEmpty()){
		   int mayor = listaFinal.get(0).size();
		   int pos =0;
		   for(int i=1; i< listaFinal.size(); i++ )
			  if(mayor<listaFinal.get(i).size()){
				 mayor = listaFinal.get(i).size();
				 pos = i;
			  }
		
		   return listaFinal.get(pos);
		}else{
		   ArrayList<Data> listNull = new ArrayList<Data>();	
		   return listNull;
		}	   
	}
	
	
	public ArrayList<Data> ComparedAtoms(boolean max, double error, Object[] atoms, Object[] atomsComp, String indice) {
		ArrayList<Data> lista = new ArrayList<Data>();
		boolean flag = true;
		while(flag){
			Data dataTemp = getMaximo();
			if(round(dataTemp.getValue().get(0),2)>=error )
				lista.add(dataTemp);
			else
				flag = false;
		}
		
        ArrayList<ArrayList<Data>> listaFinal = new ArrayList<ArrayList<Data>>();
        
        for(int i=0; i<lista.size(); i++){
        	ArrayList<Data> listaTemp = new ArrayList<Data>();
			listaTemp.add(lista.get(i));
			for(int j=0; j<lista.size(); j++){
				if(i!=j){
				   boolean flagTemp = true;
				   for(int k=0; k<listaTemp.size(); k++){
					  double valor1 = ((IAtom)atoms[lista.get(j).getPosX(0)]).getPoint3d().distance(((IAtom)atoms[listaTemp.get(k).getPosX(0)]).getPoint3d());
					  double valor2 = ((IAtom)atomsComp[lista.get(j).getPosY(0)]).getPoint3d().distance(((IAtom)atomsComp[listaTemp.get(k).getPosY(0)]).getPoint3d());
					
					  double index = 1-(Math.abs(valor1-valor2)/(valor1+valor2));
					  if(index < 0.85)
						  flagTemp = false; 
				   }
				   if(flagTemp)
					   listaTemp.add(lista.get(j));
			    }	   
		   }
			listaFinal.add(listaTemp);
		}
		
		if(!listaFinal.isEmpty()){
		   int mayor = listaFinal.get(0).size();
		   int pos =0;
		   for(int i=1; i< listaFinal.size(); i++ )
			  if(mayor<listaFinal.get(i).size()){
				 mayor = listaFinal.get(i).size();
				 pos = i;
			  }
		
		   return listaFinal.get(pos);
		}else{
		   ArrayList<Data> listNull = new ArrayList<Data>();	
		   return listNull;
		}	   
	}
	
	
	/*private ArrayList<Integer> eliminarRepetidos(ArrayList<Integer> lista){
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for(int i=0; i<lista.size(); i++){
			int pos = buscar(temp,lista.get(i));
			if(pos==-1)
				temp.add(lista.get(i));
		}
		return temp;
	}*/
	
	/*private int buscar(ArrayList<Integer> lista, int data){
		int pos =-1;
		for(int i=0; i<lista.size(); i++){
			if(lista.get(i)==data)
				pos = i;
		}
		return pos;
	}*/
	
	
	
	/*public ArrayList<Data> Compared(boolean max, double error) {
		ArrayList<Data> lista = new ArrayList<Data>();
		double mayor, menor;
		int posX = 0, posY = 0;
		boolean flag;
		if (max) {
			do {
				mayor = Double.MIN_VALUE;
				flag = false;
				for (int i = 0; i < file; i++) 
					for (int j = 0; j < column; j++) {
						@SuppressWarnings("unused")
						double valor = matrix[i][j].getValue();
						if (!matrix[i][j].getFlag())
							if ((matrix[i][j].getValue() > mayor) && (matrix[i][j].getValue() >= error)) {
								mayor = matrix[i][j].getValue();
								posX = i;
								posY = j;
								flag = true;
							}
					}
				if(flag){
					for (int i = 0; i < column; i++)
						matrix[posX][i].setFlag(true);
					for (int i = 0; i < file; i++)
						matrix[i][posY].setFlag(true);
					lista.add(matrix[posX][posY]);
				}
				
			} while (mayor >= error);
		} else {
			do {
				menor = Double.MAX_VALUE;
				flag= false;
				for (int i = 0; i < file; i++)
					for (int j = 0; j < column; j++)
						if (!matrix[i][j].getFlag()) {
							if ((matrix[i][j].getValue() < menor) && (matrix[i][j].getValue() <= error)) {
								menor = matrix[i][j].getValue();
								posX = i;
								posY = j;
								flag = true;
							}
						}
				if(flag){
					for (int i = 0; i < column; i++)
						matrix[posX][i].setFlag(true);
					for (int i = 0; i < file; i++)
						matrix[i][posY].setFlag(true);
					lista.add(matrix[posX][posY]);
				}
				
			} while (menor <= error);
		}
		return lista;
	}*/
	
	public ArrayList<Integer> ComparedDist(boolean max, double error) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		if (max) {
			ArrayList<Integer> listOut = new ArrayList<Integer>();
			/*for (int i = 0; i < file; i++) 
				 for (int j = 0; j < column; j++) 
			     if ((matrix[i][j].getValue() >= error) || i==j)
			    	 matrix[i][j].setFlag(true);
			      else
			    	 matrix[i][j].setFlag(false);*/
			
			int mayor;
			do{
				mayor= 0;
				ArrayList<Integer> fragm1 = new ArrayList<Integer>();
				@SuppressWarnings("unused")
				ArrayList<Integer> fragm2 = new ArrayList<Integer>();
				
				for (int i = 0; i < file; i++){
					int contF =0; 
					/*for (int j = 0; j < column; j++){
						if(!matrix[i][j].getFlag())  
							contF++;
					}*/
					fragm1.add(contF);
				}
				int pos =0;
				for(int i=0; i<fragm1.size(); i++)
					if(mayor < fragm1.get(i)){
						mayor = fragm1.get(i);
						pos = i;
					}
				if(mayor > 0){
					/*for(int i=0; i< file; i++){
						matrix[pos][i].setFlag(true);
						matrix[i][pos].setFlag(true);
					}*/
					listOut.add(pos);
				}					
			}while(mayor > 0);
		    for(int i=0; i< file; i++){
		    	boolean flag = false;
		    	for(int j=0; j<listOut.size(); j++)
		    		if(i == listOut.get(j))
		    			flag = true;
		    	if(!flag)
		    		list.add(i);
		    }
		} else {
			/*for (int i = 0; i < file; i++) 
			   for (int j = 0; j < column; j++) 
				  if ((matrix[i][j].getValue() <= error))
				     matrix[i][j].setFlag(true);
				  else
				     matrix[i][j].setFlag(false);*/
			
		}
		return list;
	}
	
	public ArrayList<Data> ComparedDistance(double error){
		ArrayList<Data> list = new ArrayList<Data>();
		/*for (int i = 0; i < file-1; i++)
			for (int j = i+1; j < column; j++) {
				if(matrix[i][j].getValue() < error)
					matrix[i][j].setFlag(true);
				else 
					matrix[i][j].setFlag(false);
			}*/
		
		return list;
	}
	
	public ArrayList<Data> ComparedMatix(boolean max, double error) {
		ArrayList<Data> lista = new ArrayList<Data>();
		if(max){
			for(int i=0; i < file; i++)
				for(int j=i+1; j < column; j++)
					if(matrix[i][j].getValue().get(1) >= error)
						lista.add(matrix[i][j]);
		} else {
			for(int i=0; i < file; i++)
				for(int j=i+1; j < column; j++)
					if(matrix[i][j].getValue().get(1) <= error)
						lista.add(matrix[i][j]);
		}
		return lista;
	}
	
	/*public MatrixCompared Similary(TypeSimilaryFunction funcSim, MatrixCompared A){
		SimilaryFunctions similaryFunction = new SimilaryFunctions();;
		for(int i=0; i < file; i++)
			for(int j=0; j < column; j++){
				if(i==j)
					this.addData(new Data(0.0,i,j));
				else {
					try {
						this.addData(new Data(similaryFunction.CalculeSimilaryFunction(funcSim, this.getData()[i][j].getValue(), A.getData()[i][j].getValue()),i,j));
					} catch (CDKException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
				
		return this;
	}*/
	
	public float getSimilitudPromedio(int iteraciones,boolean max) {
		int cont=0;
		double menor,mayor;
		
		if(max){
			do {
				cont++;
				mayor = Double.MIN_VALUE;
				for (int i = 0; i < file; i++) {
					for (int j = 0; j < column; j++) {
						/*if (!matrix[i][j].getFlag()) {
							valor = matrix[i][j].getValue();
							if(valor<0){
								valor *=-1;
							}
							if (valor > mayor)  {
								mayor = valor;
								posX = i;
								posY = j;
							}
						}*/
					}
				}
			  similitudPromedio +=mayor;
				
					for (int i = 0; i < column; i++) {
						//matrix[posX][i].setFlag(true);
					}				

					for (int i = 0; i < file; i++) {
						//matrix[i][posY].setFlag(true);
					}
				
			} while (cont < iteraciones);
		}
		else{
			do {
				cont++;
				menor = Double.MAX_VALUE;
				for (int i = 0; i < file; i++) {
					for (int j = 0; j < column; j++) {
						/*if (!matrix[i][j].getFlag()) {
							valor = matrix[i][j].getValue();
							if(valor<0){
								valor *=-1;
							}
							if (valor < menor)  {
								menor = valor;
								posX = i;
								posY = j;
							}
						}*/
					}
				}
			  similitudPromedio +=menor;
				
					for (int i = 0; i < column; i++) {
						//matrix[posX][i].setFlag(true);
					}				

					for (int i = 0; i < file; i++) {
						//matrix[i][posY].setFlag(true);
					}
				
			} while (cont < iteraciones);			
		}
		return (float) similitudPromedio/iteraciones;
	}
}
