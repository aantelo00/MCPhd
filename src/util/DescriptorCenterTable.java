package util;

import java.awt.event.MouseEvent;
import java.io.File;

import javax.swing.JTable;

import util.MainInterface;

public class DescriptorCenterTable extends Table {
	
	private static final long serialVersionUID = 1L;

	public DescriptorCenterTable(Object[][] data, Object[] columns, final MainInterface component, final File file){
		super(data, columns, component, file);
		setEventoMouseClicked(this);
	}		
	
	private void setEventoMouseClicked(JTable table) {
		
		table.addMouseListener(new java.awt.event.MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				molecule = buscarMoleculeActive(file);
				String cadena = getValueAt(getSelectedRow(), 1).toString();
		 		
		 		String[] temp = cadena.split(",");
		 		
		 		for (int j = 0; j < molecule.getCantAtmos(); j++) {
		 			molecule.getGraph().getAtom(j).setFlag(0, false);
		 		}
		 		
		 		
		 		for (int i = 0; i < temp.length; i++) {
					for (int j = 0; j < molecule.getCantAtmos(); j++) {
						String cadenaAtomo = molecule.getGraph().getAtom(j).getSymbol() + molecule.getGraph().getAtom(j).getID();
						if(temp[i].equals(cadenaAtomo))
							molecule.getGraph().getAtom(j).setFlag(0, true);
					}
				}
		 		
		 		//JmolPanel j = (JmolPanel)window.getToolWindowTabs()[pos].getComponent();
		 		//j.getViewer().evalString(a);
		 		//j.getViewer().evalString("set selectionhalos");	
		 		//setEnable ();
			 }
			});
	}
	
	@SuppressWarnings("unused")
	private void setEnable (){	
		  //for(int i=0; i< 4; i++)
		    //component.getRibbon().getTask(1).getBands().get(2).getControlPanel().getComponent(i).setEnabled(false);

	}
}
