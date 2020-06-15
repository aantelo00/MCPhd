package gui;

import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JLabel;

import java.awt.Font;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.JComboBox;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.DefaultComboBoxModel;
import javax.swing.filechooser.FileNameExtensionFilter;

import domain.TypeSimilaryFunction;
import util.MainInterface;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import listener.ComparedMoleculeCollectionActionListener;

public class SimilitudMCPhd {

	private JFrame frmSimilitudUtilizandoMcs;
	private JTextField textField;
	private JTextField textField_1;
	private JFileChooser jdirectorio, jdirectorio1, jdirectorio2;
	private String directorio = "";
	//private File[] files = null;
	private String directorioSalva = "";
	@SuppressWarnings("rawtypes")
	private JComboBox comboBox, comboBox_1;
	//private JTable table_1;
	private JLabel lblNewLabel, lblNewLabel_1, lblNewLabel_2, lblMtodoMcs, lblUnbral, lblActivityDirectory;
	private JTextField textField_3;
	private String index = "E";
	private JButton button, btnNewButton, btnNewButton_1, btnNewButton_2;
	private static MainInterface component = null;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					SimilitudMCPhd window = new SimilitudMCPhd();
					window.frmSimilitudUtilizandoMcs.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public SimilitudMCPhd() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void initialize() {
		frmSimilitudUtilizandoMcs = new JFrame();
		frmSimilitudUtilizandoMcs.setTitle("Similarity using MCPhd");
		frmSimilitudUtilizandoMcs.setBounds(100, 100, 613, 398);
		frmSimilitudUtilizandoMcs.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmSimilitudUtilizandoMcs.getContentPane().setLayout(null);
		
		lblNewLabel = new JLabel("Directory of mol files");
		lblNewLabel.setFont(new Font("Tahoma", Font.PLAIN, 14));
		lblNewLabel.setBounds(48, 29, 191, 16);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblNewLabel);
		
		lblNewLabel_1 = new JLabel("Output directory");
		lblNewLabel_1.setFont(new Font("Tahoma", Font.PLAIN, 14));
		lblNewLabel_1.setBounds(48, 88, 170, 16);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblNewLabel_1);
		
		lblNewLabel_2 = new JLabel("");
		lblNewLabel_2.setFont(new Font("Tahoma", Font.PLAIN, 18));
		lblNewLabel_2.setBounds(206, 294, 200, 27);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblNewLabel_2);
		
		lblMtodoMcs = new JLabel("Similarity function");
		lblMtodoMcs.setFont(new Font("Tahoma", Font.PLAIN, 14));
		lblMtodoMcs.setBounds(48, 215, 144, 16);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblMtodoMcs);
		
		lblUnbral = new JLabel("Similarity Threshold");
		lblUnbral.setFont(new Font("Tahoma", Font.PLAIN, 14));
		lblUnbral.setBounds(262, 215, 144, 16);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblUnbral);
		
		lblActivityDirectory = new JLabel("Activity file");
		lblActivityDirectory.setFont(new Font("Tahoma", Font.PLAIN, 14));
		lblActivityDirectory.setBounds(48, 149, 121, 16);
		frmSimilitudUtilizandoMcs.getContentPane().add(lblActivityDirectory);
		
		textField = new JTextField();
		textField.setEditable(false);
		textField.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textField.setBounds(48, 51, 435, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(textField);
		textField.setColumns(10);
		
		textField_1 = new JTextField();
		textField_1.setEditable(false);
		textField_1.setEnabled(false);
		textField_1.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textField_1.setBounds(48, 112, 435, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(textField_1);
		textField_1.setColumns(10);
		
		textField_3 = new JTextField();
		textField_3.setEditable(false);
		textField_3.setFont(new Font("Tahoma", Font.PLAIN, 14));
		textField_3.setEnabled(false);
		textField_3.setBounds(48, 172, 435, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(textField_3);
		textField_3.setColumns(10);
		
		comboBox_1 = new JComboBox();
		comboBox_1.setEnabled(false);
		comboBox_1.setModel(new DefaultComboBoxModel(TypeSimilaryFunction.values()));
		comboBox_1.setFont(new Font("Tahoma", Font.PLAIN, 14));
		comboBox_1.setBounds(48, 237, 191, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(comboBox_1);
		
		comboBox = new JComboBox();
		comboBox.setEnabled(false);
		comboBox.setModel(new DefaultComboBoxModel(new String[] {"0.95"}));
		comboBox.setFont(new Font("Tahoma", Font.PLAIN, 14));
		comboBox.setBounds(262, 237, 119, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(comboBox);
		
		button = new JButton("...");
		button.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				
				FileNameExtensionFilter xmlFilter = new FileNameExtensionFilter("mol files (*.mol)", "mol");
				jdirectorio = new JFileChooser();
	    		jdirectorio.setApproveButtonText("Seleccionar Directorio");
	    		jdirectorio.setDialogTitle("Seleccionar el directorio.");
	    		jdirectorio.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	    		jdirectorio.addChoosableFileFilter(xmlFilter);
	    		jdirectorio.setFileFilter(xmlFilter);
	    		jdirectorio.showOpenDialog(new JFrame());
	    		if(jdirectorio.getSelectedFile() != null){
	    			directorio = jdirectorio.getSelectedFile().getPath();
	    		    textField.setText(directorio);
	    		    textField_1.setEnabled(true);
	    		    btnNewButton.setEnabled(true);
	    		    lblNewLabel_2.setText("");
	    		}    
			}
		});
		button.setFont(new Font("Tahoma", Font.PLAIN, 14));
		button.setBounds(498, 50, 58, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(button);
		
		btnNewButton = new JButton("...");
		btnNewButton.setEnabled(false);
		btnNewButton.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				jdirectorio1 = new JFileChooser();
	    		jdirectorio1.setApproveButtonText("Seleccionar Directorio");
	    		jdirectorio1.setDialogTitle("Seleccionar el directorio.");
	    		jdirectorio1.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	    		jdirectorio1.setSelectedFile(null);
	    		jdirectorio1.showOpenDialog(new JFrame());
	    		if(jdirectorio1.getSelectedFile() != null){
	    			directorioSalva = jdirectorio1.getSelectedFile().getPath();
	    		    textField_1.setText(directorioSalva);
	    		    textField_3.setEnabled(true);
	    		    btnNewButton_2.setEnabled(true);
	    		}
			}
		});
		btnNewButton.setFont(new Font("Tahoma", Font.PLAIN, 14));
		btnNewButton.setBounds(498, 111, 58, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(btnNewButton);
		
		btnNewButton_1 = new JButton("Calculate similarity");
		btnNewButton_1.addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent arg0) {
				
				TypeSimilaryFunction simFunction = (TypeSimilaryFunction)comboBox_1.getSelectedItem();
				double error = Double.parseDouble(comboBox.getSelectedItem().toString());
				ComparedMoleculeCollectionActionListener cal = new ComparedMoleculeCollectionActionListener(component, textField.getText(), textField_1.getText(), textField_3.getText(), simFunction , error, index);
			    try {
					cal.comparedMolecular();
					lblNewLabel_2.setText("end of processing");
					btnNewButton_1.setEnabled(false);
					textField.setText("");
					textField_1.setText("");
					textField_3.setText("");
					comboBox.setEnabled(false);
					comboBox_1.setEnabled(false);
					btnNewButton.setEnabled(false);
					btnNewButton_2.setEnabled(false);
					directorio = "";
				} catch (IOException e) {
					Object[] options = { "Ok" };
					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);
				} catch (CloneNotSupportedException e) {
					Object[] options = { "Ok" };
					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);
				} catch (Exception e){
					Object[] options = { "Ok" };
					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);
				}
			}
		});
		btnNewButton_1.setEnabled(false);
		btnNewButton_1.setFont(new Font("Tahoma", Font.PLAIN, 14));
		btnNewButton_1.setBounds(412, 237, 144, 24);
		frmSimilitudUtilizandoMcs.getContentPane().add(btnNewButton_1);
		
		btnNewButton_2 = new JButton("...");
		btnNewButton_2.setEnabled(false);
		btnNewButton_2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				FileNameExtensionFilter xmlFilter = new FileNameExtensionFilter("csv files (*.csv)", "csv");


				jdirectorio2 = new JFileChooser();
	    		jdirectorio2.setApproveButtonText("Seleccionar Directorio");
	    		jdirectorio2.setDialogTitle("Seleccionar el directorio.");
	    		jdirectorio2.setFileSelectionMode(JFileChooser.FILES_ONLY);
	    		jdirectorio2.addChoosableFileFilter(xmlFilter);
	    		jdirectorio2.setFileFilter(xmlFilter);
	    		jdirectorio2.showOpenDialog(new JFrame()); 		
	    		if(jdirectorio2.getSelectedFile() != null){
	    			directorio = jdirectorio2.getSelectedFile().getPath();
	    		    textField_3.setText(directorio);
	    		    comboBox_1.setEnabled(true);
	    		    comboBox.setEnabled(true);
	    		    btnNewButton_1.setEnabled(true);
		         }   
			}
		});
		btnNewButton_2.setBounds(498, 171, 58, 25);
		frmSimilitudUtilizandoMcs.getContentPane().add(btnNewButton_2);
	}
}
