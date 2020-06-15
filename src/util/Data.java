package util;

import java.util.ArrayList;

public class Data {
	
	private ArrayList<Double> value;
	private int valor;
	private ArrayList<Integer> posX = new ArrayList<Integer>();;
	private ArrayList<Integer> posY = new ArrayList<Integer>();
	private boolean flag;
	
	public Data(ArrayList<Double> value, ArrayList<Integer> posX, ArrayList<Integer> posY){
		this.value = value;
		this.posX = posX;
		this.posY = posY;
		this.flag = false;
	}
	
	public Data(int valor){
		this.valor = valor;
	}
	
	public ArrayList<Double> getValue(){
		return value;
	}
	
	public int getValor(){
		return valor;
	}
	
	public int getPosX(int pos){
		return posX.get(pos);
	}
	
	public int getPosY(int pos){
		return posY.get(pos);
	}
	
	public boolean getFlag(){
		return flag;
	}
	
	public void setFlag(boolean valor){
		this.flag = valor;
	}
	
}
