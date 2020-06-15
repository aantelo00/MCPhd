package util;

import javax.swing.*;

import domain.TypeSimilaryFunction;

public class SearchOptions {
	public TypeSimilaryFunction simFunction;
	public double error;
	public boolean isMaximal;	

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public SearchOptions() {
		JPanel myPanel = new JPanel();
		JLabel simFuncLabel = new JLabel("Funcion:");
		JComboBox simFunc = new JComboBox();
		JLabel errorlabel = new JLabel("Error:");
		JTextField error = new JTextField(10);
		JCheckBox max = new JCheckBox("Es máximo?");

		simFunc.setModel(new DefaultComboBoxModel(TypeSimilaryFunction.values()));

		myPanel.add(simFuncLabel);
		myPanel.add(simFunc);
		myPanel.add(Box.createHorizontalStrut(15));
		myPanel.add(errorlabel);
		myPanel.add(error);
		myPanel.add(Box.createHorizontalStrut(15));
		myPanel.add(max);		

		int result = JOptionPane.showConfirmDialog(null, myPanel, "Introduzca las opciones de búsqueda", JOptionPane.OK_CANCEL_OPTION);
		if (result == JOptionPane.OK_OPTION) {
			this.simFunction = (TypeSimilaryFunction)simFunc.getSelectedItem();
			this.isMaximal = max.isSelected();
			this.error = Double.parseDouble(error.getText());
		}
	}

}
