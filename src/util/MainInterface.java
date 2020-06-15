package util;

import java.util.ArrayList;

import javax.swing.JTree;

//import util.Normalization;

import domain.Molecule;


public interface MainInterface{
	public ArrayList<String> getMySources();
	public ArrayList<String> getMySourcesPath();
	public ArrayList<String> getRecientes();
	public JTree getTree();
	public void UpdateTree();
	public ArrayList<String> getMySourcesCompared();
	public ArrayList<Molecule> getListMolecule();
	//public ArrayList<Normalization> getNormalization();
	public boolean getShowHidrogens();
	public boolean getShowHalos();
	public void setShowHidrogens(boolean value);
	public void setShowHalos(boolean value);
}
