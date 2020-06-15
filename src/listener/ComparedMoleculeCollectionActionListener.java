package listener;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.jmol.JmolPanel;
import org.openscience.cdk.exception.CDKException;

import util.MainInterface;
import domain.FragmentoMolecular;
import domain.FragmentosPMC;
import domain.Molecule;
import domain.ReducedGraph;
import domain.TypeSimilaryFunction;

public class ComparedMoleculeCollectionActionListener {
	
	private MainInterface component = null;
	private String index;
	private String entryDirectory, outputDirectory, activityFile;
	TypeSimilaryFunction simFunction;
	double threshold;
	@SuppressWarnings("unused")
	private Molecule molecule, moleculeTemp, moleculeCompared;
	private JFileChooser jdirectorio;
	private NumberFormat format = NumberFormat.getInstance();
	private String directorio = "";
	
	public ComparedMoleculeCollectionActionListener(MainInterface component, String entryDirectory, String outputDirectory, String activityFile, TypeSimilaryFunction simFunction, double threshold, String index) {
		this.component = (MainInterface) component;
		this.index = index;
		this.entryDirectory = entryDirectory;
		this.simFunction = simFunction;
		this.threshold = threshold;
		this.outputDirectory = outputDirectory;
		this.activityFile = activityFile;
		jdirectorio = new JFileChooser();
		format.setMaximumFractionDigits(5);
	}
	
public void comparedMolecular() throws IOException, CloneNotSupportedException {
		
	    
		//molecule = buscarMoleculeActive();
        
       // if (molecule != null) {
        	
        	File[] files = null;
   
        	 FilenameFilter filterMOL = new FilenameFilter() {
     			public boolean accept(File dir, String name) {
     				return name.endsWith(".mol");
     			}
     		};
     		
     	    files = getMoleculas(new File(entryDirectory), filterMOL);
     	    String moleculeTempActivity="", moleculeComparedActivity="";
    		for(int i=0; i< files.length; i++){
    			moleculeTemp = new Molecule();
    			try {
    				moleculeTemp.setMolecula(files[i]);
    				String FicheroSalva = outputDirectory + System.getProperty("file.separator") + moleculeTemp.getNombre().substring(0, moleculeTemp.getNombre().length()-4) +".csv";
    				FileWriter salva = new FileWriter(FicheroSalva);
    				salva.write("Molecule 1;Activity;Fragment;Molecule 2;Activity;Fragment;Similarity" + "\n");
					moleculeTemp.calcularIndicesTopograficos();
					ReducedGraph graph = moleculeTemp.getReducedGraph();
					moleculeTemp.calcularDistTopograficas(graph); 
					try {
						moleculeTempActivity = searchActivity(moleculeTemp.getNombre());
					} catch (Exception e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					} 

					
	    			for (int k=0; k< files.length; k++){
	    				System.out.println(k + "- " + files[k].getName());
	    				String cadena = "";
	    				FragmentoMolecular f1 = new FragmentoMolecular();
	    				FragmentoMolecular f2 = new FragmentoMolecular();
	        			moleculeCompared = new Molecule();
	        			try {
	        				moleculeCompared.setMolecula(files[k]);
	        				moleculeCompared.calcularIndicesTopograficos();
	        				ReducedGraph graphTemp = moleculeCompared.getReducedGraph();
	    					moleculeCompared.calcularDistTopograficas(graphTemp);
	        			} catch (Exception e) {
	        				Object[] options = { "Ok" };
	    					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);
	        				//System.out.println(files[k].getName());
	        			}	
	        			moleculeTemp.normalizationDescriptor(moleculeCompared);
	        			
	        			try {
							moleculeComparedActivity = searchActivity(moleculeCompared.getNombre());
						} catch (Exception e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
	        			
	        			FragmentosPMC fragmento;
	        			
	        			@SuppressWarnings("unused")
	        			ArrayList<FragmentosPMC> resultados = new ArrayList<FragmentosPMC>();
	        			int maxGrado = 0;
	        			@SuppressWarnings("unused")
						JmolPanel jmp = new JmolPanel();
	        			double indexTemp = 0.0;
	        			try {
	        				fragmento = moleculeTemp.getFragmentoPMCCD(simFunction, threshold, moleculeCompared, index);
	        				if (fragmento != null) {
	        					int grado = fragmento.getFragMol1().getGradoFragmento();
	        					if (grado > maxGrado)
	        						maxGrado = grado;
	        					f1 = fragmento.getFragMol1();
	        					f2 = fragmento.getFragMol2();
	        					
	        					//double MCPA = f1.getValoresTopograficosPMCCD(index);
	        					//double MCPB = f2.getValoresTopograficosPMCCD(index);
	        					//double MCP = Math.min(MCPA, MCPB);
	        					//double TPA = moleculeTemp.sumaTotalDescriptorTopographicalNorm(index);
	        					//double TPB = moleculeCompared.sumaTotalDescriptorTopographicalNorm(index);
	        					
	        					//double c = (f1.getListAtomsCD().size()+f2.getListAtomsCD().size())/2;
	        					double c = Math.min(f1.getListAtomsCD().size(),f2.getListAtomsCD().size());
	        					double a = moleculeTemp.getCantAtmosPesados();
	        					double b = moleculeCompared.getCantAtmosPesados();
	        					
	        					switch (simFunction) {
	     					         case TANIMOTO: {
	     						            indexTemp = c/(a+b-c);
	     					        	    //indexTemp = MCP/(TPA+TPB-MCP);
	     						            break;
	     					         }
	     					         /*case DICE: {
	     						            //indexTemp = (PMC)/(0.5*(PTA)+(PTB));
	     						            break;
	     					         }
	     					         case COSINE: {
	     						            //indexTemp = (PMC)/(Math.sqrt((PTA)*(PTB)));
	     						            break;
	     					          }*/
	     					     }
	        					cadena += format.format(indexTemp);
	     						salva.write(moleculeTemp.getNombre() + ";" + moleculeTempActivity + ";" + f1.getTipo() + ";" + moleculeCompared.getNombre() + ";" + moleculeComparedActivity 
	            					    + ";" + f2.getTipo() + ";" +cadena + "\n");
	        				}
	        				else{
	        					cadena = "0,00000;";
	        					salva.write(moleculeTemp.getNombre() + ";" + moleculeTempActivity + ";" + "" + ";" + moleculeCompared.getNombre() + ";" + moleculeComparedActivity 
	            					    + ";" + "" + ";" + cadena + "\n");	
	        				}
	        			} catch (CDKException e) {
	        				Object[] options = { "Ok" };
	    					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);

	        			}
	        			
	        			
	        		}
	    			System.out.println(moleculeTemp.getNombre());
	    			salva.close();
				} catch (CDKException e1) {
					Object[] options = { "Ok" };
					JOptionPane.showOptionDialog(new JFrame(), "I can't perform the calculations, check the directory of the molecules", "Error!!", JOptionPane.YES_OPTION, JOptionPane.ERROR_MESSAGE, null, options, options[0]);

				}
    			
           }
    	    System.out.println("Finish");
}

private  File[] getMoleculas(File dir, FilenameFilter filter){
	File[] result = null;
	result = dir.listFiles(filter);	
	return result;
}

private String searchActivity(String name) throws Exception, IOException, FileNotFoundException{
	String activity = "";
	FileReader fileAct;
	try {
		fileAct = new FileReader(activityFile);
		BufferedReader readAct = new BufferedReader(fileAct);
		boolean inicio = true;
		//readAct.readLine();
		String cadena = "";
		while (inicio) {
			cadena = readAct.readLine();
			String[] cadenaTemp = cadena.split(";");
			if (cadenaTemp[1].equals(name.substring(0, name.length()-4))) {
				activity = cadenaTemp[2];
				inicio = false;
				readAct.close();
				fileAct.close();
			}
		}
	} catch (FileNotFoundException e) {
		Logger.getLogger(ComparedMoleculeCollectionActionListener.class.getName()).log(Level.SEVERE, null, e);
	} catch (IOException e) {
		Logger.getLogger(ComparedMoleculeCollectionActionListener.class.getName()).log(Level.SEVERE, null, e);
	} catch (Exception e){
		Logger.getLogger(ComparedMoleculeCollectionActionListener.class.getName()).log(Level.SEVERE, null, e);	}
	return activity;
}

}
