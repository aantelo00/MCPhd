package util;

import java.io.File;

import javax.swing.JTable;

import util.MainInterface;

import domain.Molecule;

@SuppressWarnings("serial")
public abstract class Table extends JTable {
	protected File file;
	protected MainInterface component;
	protected Molecule molecule;
	
	public Table(Object[][] data, Object[] columns, MainInterface component, final File file){
		super(data, columns);
		this.file = file;
		this.component = component;		
		
	}	
	
	public Molecule buscarMoleculeActive(File file) {
		Molecule molecule = null;
		for (int i = 0; i < component.getListMolecule().size(); i++) {
			if(file.getName().equals(component.getListMolecule().get(i).getGraph().getID())){
				molecule = component.getListMolecule().get(i);
			}
				
		}
		
		return molecule;
	}

}
