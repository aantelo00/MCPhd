package util;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.DataFormat;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFRichTextString;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;

import domain.CentroDescriptor;
import domain.FragmentoMolecular;
import domain.FragmentosPMC;
import domain.TypeSimilaryFunction;

public class ExcelStadistics {
	private XSSFWorkbook book;
	private XSSFSheet sheet;

	public ExcelStadistics() {
		this.book = new XSSFWorkbook();
		this.sheet = book.createSheet("Resultados");
	}

	public void createStadisticsFragments(ArrayList<FragmentoMolecular> fragments, String path, char indice, TypeSimilaryFunction funcSim, double error) throws CDKException {
		String[] headerFragment = null;
		switch(indice){
		    case 'A': headerFragment = new String[] { "GR�FICO", "FRAGMENTO", "MOL�CULA", "ETPG-CD1", "DIF", "RTPG-CD1", "DIF", "LTPG-CD1", "DIF", "ETPG-CD2", "DIF", "RTPG-CD2", "DIF", "LTPG-CD2", "DIF", "DISTANCIA", "DIF", "ACTIVIDAD", "ENSAYO", "TANIMOTO"};
                      break;
            case 'E': headerFragment = new String[] { "GR�FICO", "FRAGMENTO", "MOL�CULA", "ETPG-CD1", "DIF", "ETPG-CD2", "DIF", "DISTANCIA", "DIF", "ACTIVIDAD", "ENSAYO", "TANIMOTO"};
                      break;
            case 'R': headerFragment = new String[] { "GR�FICO", "FRAGMENTO", "MOL�CULA", "RTPG-CD1", "DIF", "RTPG-CD2", "DIF", "DISTANCIA", "DIF", "ACTIVIDAD", "ENSAYO", "TANIMOTO"};
                      break;
            case 'L': headerFragment = new String[] { "GR�FICO", "FRAGMENTO", "MOL�CULA", "LTPG-CD1", "DIF", "LTPG-CD2", "DIF", "DISTANCIA", "DIF", "ACTIVIDAD", "ENSAYO", "TANIMOTO"};
                      break;
		}
		
		createHeader(headerFragment);
		int row = 2;
		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];
		FragmentoMolecular target = fragments.get(0);
		for (FragmentoMolecular f : fragments) {
			Row bodyRow = sheet.createRow(row++);
			CreateCell(bodyRow, bodyStyle(), 0, "");
			CreateCell(bodyRow, bodyStyle(), 1, f.getTipo());
			CreateCell(bodyRow, bodyStyle(), 2, f.getMolecula().getNombre());
			
			switch(indice){
			   case 'A': {
				    CreateCell(bodyRow, bodyStyle(), 3, f.vectorValoresTopograficos()[0]);
				    CreateCell(bodyRow, bodyStyle(), 4, valueAbs(target.vectorValoresTopograficos()[0], f.vectorValoresTopograficos()[0]));
				    CreateCell(bodyRow, bodyStyle(), 5, f.vectorValoresTopograficos()[1]);
					CreateCell(bodyRow, bodyStyle(), 6, valueAbs(target.vectorValoresTopograficos()[1], f.vectorValoresTopograficos()[1]));
					CreateCell(bodyRow, bodyStyle(), 7, f.vectorValoresTopograficos()[2]);
					CreateCell(bodyRow, bodyStyle(), 8, valueAbs(target.vectorValoresTopograficos()[2], f.vectorValoresTopograficos()[2]));
					CreateCell(bodyRow, bodyStyle(), 9, f.vectorValoresTopograficos()[3]);
					CreateCell(bodyRow, bodyStyle(), 10, valueAbs(target.vectorValoresTopograficos()[3], f.vectorValoresTopograficos()[3]));
					CreateCell(bodyRow, bodyStyle(), 11, f.vectorValoresTopograficos()[4]);
					CreateCell(bodyRow, bodyStyle(), 12, valueAbs(target.vectorValoresTopograficos()[4], f.vectorValoresTopograficos()[4]));
					CreateCell(bodyRow, bodyStyle(), 13, f.vectorValoresTopograficos()[5]);
					CreateCell(bodyRow, bodyStyle(), 14, valueAbs(target.vectorValoresTopograficos()[5], f.vectorValoresTopograficos()[5]));
					CreateCell(bodyRow, bodyStyle(), 15, f.vectorValoresTopograficos()[6]);
					CreateCell(bodyRow, bodyStyle(), 16, valueAbs(target.vectorValoresTopograficos()[6], f.vectorValoresTopograficos()[6]));
					CreateCell(bodyRow, bodyStyle(), 17, f.getMolecula().getActividad());
					CreateCell(bodyRow, bodyStyle(), 18, f.getMolecula().getDireccion().split(File.separator + File.separator)[f.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);
					CreateCell(bodyRow, bodyStyle(), 19, f.getIndexSimilary());
					break;
			   }
			   case 'E': {
				    CreateCell(bodyRow, bodyStyle(), 3, f.vectorValoresTopograficos()[0]);
				    CreateCell(bodyRow, bodyStyle(), 4, valueAbs(target.vectorValoresTopograficos()[0], f.vectorValoresTopograficos()[0]));
				    CreateCell(bodyRow, bodyStyle(), 5, f.vectorValoresTopograficos()[3]);
				    CreateCell(bodyRow, bodyStyle(), 6, valueAbs(target.vectorValoresTopograficos()[3], f.vectorValoresTopograficos()[3]));
					CreateCell(bodyRow, bodyStyle(), 7, f.vectorValoresTopograficos()[6]);
					CreateCell(bodyRow, bodyStyle(), 8, valueAbs(target.vectorValoresTopograficos()[6], f.vectorValoresTopograficos()[6]));
					CreateCell(bodyRow, bodyStyle(), 9, f.getMolecula().getActividad());
					CreateCell(bodyRow, bodyStyle(), 10, f.getMolecula().getDireccion().split(File.separator + File.separator)[f.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);
					CreateCell(bodyRow, bodyStyle(), 11, f.getIndexSimilary());
					break;
			   }
			   case 'R': {
				    CreateCell(bodyRow, bodyStyle(), 3, f.vectorValoresTopograficos()[1]);
				    CreateCell(bodyRow, bodyStyle(), 4, valueAbs(target.vectorValoresTopograficos()[1], f.vectorValoresTopograficos()[1]));
				    CreateCell(bodyRow, bodyStyle(), 5, f.vectorValoresTopograficos()[4]);
				    CreateCell(bodyRow, bodyStyle(), 6, valueAbs(target.vectorValoresTopograficos()[4], f.vectorValoresTopograficos()[4]));
					CreateCell(bodyRow, bodyStyle(), 7, f.vectorValoresTopograficos()[6]);
					CreateCell(bodyRow, bodyStyle(), 8, valueAbs(target.vectorValoresTopograficos()[6], f.vectorValoresTopograficos()[6]));
					CreateCell(bodyRow, bodyStyle(), 9, f.getMolecula().getActividad());
					CreateCell(bodyRow, bodyStyle(), 10, f.getMolecula().getDireccion().split(File.separator + File.separator)[f.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);
					CreateCell(bodyRow, bodyStyle(), 11, f.getIndexSimilary());
					break;
			   }
			   case 'L': {
				    CreateCell(bodyRow, bodyStyle(), 3, f.vectorValoresTopograficos()[2]);
				    CreateCell(bodyRow, bodyStyle(), 4, valueAbs(target.vectorValoresTopograficos()[2], f.vectorValoresTopograficos()[2]));
				    CreateCell(bodyRow, bodyStyle(), 5, f.vectorValoresTopograficos()[5]);
				    CreateCell(bodyRow, bodyStyle(), 6, valueAbs(target.vectorValoresTopograficos()[5], f.vectorValoresTopograficos()[5]));
					CreateCell(bodyRow, bodyStyle(), 7, f.vectorValoresTopograficos()[6]);
					CreateCell(bodyRow, bodyStyle(), 8, valueAbs(target.vectorValoresTopograficos()[6], f.vectorValoresTopograficos()[6]));
					CreateCell(bodyRow, bodyStyle(), 9, f.getMolecula().getActividad());
					CreateCell(bodyRow, bodyStyle(), 10, f.getMolecula().getDireccion().split(File.separator + File.separator)[f.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);
					CreateCell(bodyRow, bodyStyle(), 11, f.getIndexSimilary());
					break;
			   }
			}
		}
		saveDocument(collection, target.getTipo(), funcSim.toString(), error);
	}
	
	public void createMatrixSimilariry(ArrayList<IAtom> listAtoms, ArrayList<IAtom> listAtomsComp, MatrixCompared matrix, String path) throws CDKException {
		String[] headerFragment = new String[listAtomsComp.size()+1];
		headerFragment[0]="";
		for (int i=0; i < listAtomsComp.size(); i++) {
			headerFragment[i+1] = listAtomsComp.get(i).getSymbol() + listAtomsComp.get(i).getID();
		}
		createHeader(headerFragment);
		int row = 2;
		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];
		for(int i=0; i < matrix.getFile(); i++){
			Row bodyRow = sheet.createRow(row++);
		    CreateCell(bodyRow, bodyStyle(), 0, listAtoms.get(i).getSymbol() + listAtoms.get(i).getID());
			for(int j=0; j < matrix.getColumn(); j++){
				CreateCell(bodyRow, bodyStyle(), j+1, matrix.getData()[i][j].getValue().get(0));
			}
		}	
		
		saveDocument(collection, "", "", 0.0);
	}
	
	public void createMatrixSimilariry(Object[] centroDescriptores, Object[] centroDescriptoresComp, MatrixCompared matrix, String path) throws CDKException {
		String[] headerFragment = new String[centroDescriptoresComp.length+1];
		headerFragment[0]="";
		for (int i=0; i < centroDescriptoresComp.length; i++) {
			headerFragment[i+1] = (((CentroDescriptor)centroDescriptoresComp[i])).getNombreCD();
		}
		createHeader(headerFragment);
		int row = 2;
		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];
		for(int i=0; i < matrix.getFile(); i++){
			Row bodyRow = sheet.createRow(row++);
		    CreateCell(bodyRow, bodyStyle(), 0, (((CentroDescriptor)centroDescriptores[i])).getNombreCD());
			for(int j=0; j < matrix.getColumn(); j++){
				CreateCell(bodyRow, bodyStyle(), j+1, matrix.getData()[i][j].getValue().get(0));
			}
		}	
		
		saveDocument(collection, "", "", 0.0);
	}
	
	public void createMatrixSimilariryAtoms(Object[] atoms, Object[] atomsComp, MatrixCompared matrix, String path) throws CDKException {
		String[] headerFragment = new String[atomsComp.length+1];
		headerFragment[0]="";
		for (int i=0; i < atomsComp.length; i++) {
			headerFragment[i+1] = (((IAtom)atomsComp[i])).getSymbol()+ (((IAtom)atomsComp[i])).getID() ;
		}
		createHeader(headerFragment);
		int row = 2;
		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];
		for(int i=0; i < matrix.getFile(); i++){
			Row bodyRow = sheet.createRow(row++);
		    CreateCell(bodyRow, bodyStyle(), 0, (((IAtom)atoms[i])).getSymbol()+(((IAtom)atoms[i])).getID());
			for(int j=0; j < matrix.getColumn(); j++){
				CreateCell(bodyRow, bodyStyle(), j+1, matrix.getData()[i][j].getValue().get(0));
			}
		}	
		
		saveDocument(collection, "", "", 0.0);
	}
	
	public void createMatrixSimilariry(Object[] centroDescriptores, Object[] centroDescriptoresComp, Data[][] matrix, String path) throws CDKException {
		String[] headerFragment = new String[centroDescriptoresComp.length+1];
		headerFragment[0]="";
		for (int i=0; i < centroDescriptoresComp.length; i++) {
			headerFragment[i+1] = (((CentroDescriptor)centroDescriptoresComp[i])).getNombreCD();
		}
		createHeader(headerFragment);
		int row = 2;
		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];
		for(int i=0; i < matrix.length; i++){
			Row bodyRow = sheet.createRow(row++);
		    CreateCell(bodyRow, bodyStyle(), 0, (((CentroDescriptor)centroDescriptores[i])).getNombreCD());
			for(int j=0; j < matrix[0].length; j++){
				CreateCell(bodyRow, bodyStyle(), j+1, matrix[i][j].getValue().get(0));
			}
		}	
		
		saveDocument(collection, "", "", 0.0);
	}

	public void createStadisticsMolecules(ArrayList<FragmentosPMC> fragments, String path, TypeSimilaryFunction funcSim, int maxGrado, double error) {
		String[] headerMolecule = new String[maxGrado * 6 + 8];

		headerMolecule[0] = "GRAFICO";
		headerMolecule[1] = "FRAGMENTO";
		headerMolecule[2] = "MOLECULE";
		headerMolecule[3] = "ENSAYO";
		int p = 4;
		for (int i = 0; i < maxGrado; i++) {
			headerMolecule[p++] = "ETPG-CD" + (i + 1);
			headerMolecule[p++] = "DIF";
			headerMolecule[p++] = "RTPG-CD" + (i + 1);
			headerMolecule[p++] = "DIF";
			headerMolecule[p++] = "LTPG-CD" + (i + 1);
			headerMolecule[p++] = "DIF";
		}
		headerMolecule[p++] = "DISTANCIA MEDIA";
		headerMolecule[p++] = "DIF";
		headerMolecule[p++] = "ACTIVIDAD";
		headerMolecule[p++] = "SIMILITUD";

		createHeader(headerMolecule);
		int row = 2;

		String collection = path.split(File.separator + File.separator)[path.split(File.separator + File.separator).length - 1];

		for (FragmentosPMC f : fragments) {
			Row bodyRow = sheet.createRow(row++);

			FragmentoMolecular frg1 = f.getFragMol1();
			FragmentoMolecular frg2 = f.getFragMol2();

			CreateCell(bodyRow, bodyStyle(), 0, "");
			CreateCell(bodyRow, bodyStyle(), 1, frg1.getTipo());
			CreateCell(bodyRow, bodyStyle(), 2, frg1.getMolecula().getNombre());
			CreateCell(bodyRow, bodyStyle(), 3, frg1.getMolecula().getDireccion().split(File.separator + File.separator)[frg1.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);
			p = 4;
			for (int i = 0, j = 0; i < maxGrado; i++) {
				if (i < frg1.getGradoFragmento()) {
					CreateCell(bodyRow, bodyStyle(), p++, frg1.vectorValoresTopograficos('A')[j++]);
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, frg1.vectorValoresTopograficos('A')[j++]);
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, frg1.vectorValoresTopograficos('A')[j++]);
					CreateCell(bodyRow, bodyStyle(), p++, "");
				} else {
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
				}
			}
			CreateCell(bodyRow, bodyStyle(), p++, frg1.vectorValoresTopograficos('A')[frg1.getGradoFragmento() * 3]);
			CreateCell(bodyRow, bodyStyle(), p++, "");
			CreateCell(bodyRow, bodyStyle(), p++, frg1.getMolecula().getActividad());
			CreateCell(bodyRow, bodyStyle(), p++, "");

			bodyRow = sheet.createRow(row++);

			CreateCell(bodyRow, bodyStyle(), 0, "");
			CreateCell(bodyRow, bodyStyle(), 1, frg2.getTipo());
			CreateCell(bodyRow, bodyStyle(), 2, frg2.getMolecula().getNombre());
			CreateCell(bodyRow, bodyStyle(), 3, frg2.getMolecula().getDireccion().split(File.separator + File.separator)[frg2.getMolecula().getDireccion().split(File.separator + File.separator).length - 2]);

			p = 4;
			for (int i = 0, j = 0; i < maxGrado; i++) {
				if (i < frg2.getGradoFragmento()) {
					CreateCell(bodyRow, bodyStyle(), p++, frg2.vectorValoresTopograficos('A')[j]);
					CreateCell(bodyRow, bodyStyle(), p++, valueAbs(frg1.vectorValoresTopograficos('A')[j], frg2.vectorValoresTopograficos('A')[j++]));
					CreateCell(bodyRow, bodyStyle(), p++, frg2.vectorValoresTopograficos('A')[j]);
					CreateCell(bodyRow, bodyStyle(), p++, valueAbs(frg1.vectorValoresTopograficos('A')[j], frg2.vectorValoresTopograficos('A')[j++]));
					CreateCell(bodyRow, bodyStyle(), p++, frg2.vectorValoresTopograficos('A')[j]);
					CreateCell(bodyRow, bodyStyle(), p++, valueAbs(frg1.vectorValoresTopograficos('A')[j], frg2.vectorValoresTopograficos('A')[j++]));
				} else {
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
					CreateCell(bodyRow, bodyStyle(), p++, "");
				}
			}

			CreateCell(bodyRow, bodyStyle(), p++, frg2.vectorValoresTopograficos('A')[frg2.getGradoFragmento() * 3]);
			CreateCell(bodyRow, bodyStyle(), p++, valueAbs(frg1.vectorValoresTopograficos('A')[frg2.getGradoFragmento() * 3], frg2.vectorValoresTopograficos('A')[frg2.getGradoFragmento() * 3]));
			CreateCell(bodyRow, bodyStyle(), p++, frg2.getMolecula().getActividad());
			CreateCell(bodyRow, bodyStyle(), p++, "");
			row++;
		}

		saveDocument(collection, fragments.get(0).getFragMol1().getMolecula().getNombre(), funcSim.toString(), error);
	}

	
	
	
	private double valueAbs(double va1, double val2) {
		return Math.abs(Math.abs(va1) - Math.abs(val2));
	}

	private void CreateCell(Row parentRow, CellStyle cstyle, int number, String value) {
		Cell celda = parentRow.createCell((short) number);
		celda.setCellValue(value);
		celda.setCellStyle(cstyle);
	}

	private void CreateCell(Row parentRow, CellStyle cstyle, int number, double value) {
		Cell celda = parentRow.createCell((short) number);
		celda.setCellValue(value);
		celda.setCellStyle(cstyle);
	}

	private void createHeader(String[] header) {
		int columns = header.length;
		Row headerRow = sheet.createRow(1);
		headerRow.setHeightInPoints(30);

		for (int i = 0; i < columns; i++) {
			Cell celda = headerRow.createCell((short) i);
			XSSFRichTextString texto = new XSSFRichTextString(header[i]);
			celda.setCellValue(texto);
			celda.setCellStyle(headerStyle());
		}
	}

	public void saveDocument(String collection, String target, String funcSim, double error) {
		try {
			JFileChooser fc = new JFileChooser();
			fc.setAcceptAllFileFilterUsed(false);
			fc.setFileFilter(new FileNameExtensionFilter("Libro de Excel (*.xlsx)", "xlsx"));
			fc.setSelectedFile(new File(System.getProperty("user.home") + "\\Busqueda_" + funcSim + "_" + error + "_" + target + "_" + collection + ".xlsx"));
			int returnVal = fc.showSaveDialog(null);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				String dir = file.getAbsolutePath();
				if (dir.lastIndexOf(".") != -1) {
					dir = dir.substring(0, dir.lastIndexOf("."));
				}

				book.write(new FileOutputStream(new File(dir + ".xlsx")));

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public CellStyle headerStyle() {
		CellStyle cs = book.createCellStyle();
		cs.setFillForegroundColor(IndexedColors.GREY_25_PERCENT.getIndex());
		cs.setFillPattern(CellStyle.SOLID_FOREGROUND);
		cs.setAlignment(CellStyle.ALIGN_CENTER);
		cs.setVerticalAlignment(CellStyle.VERTICAL_CENTER);
		cs.setBorderLeft(CellStyle.BORDER_THIN);
		cs.setBorderTop(CellStyle.BORDER_THIN);
		cs.setBorderRight(CellStyle.BORDER_THIN);
		cs.setBorderBottom(CellStyle.BORDER_THIN);
		cs.setWrapText(true);

		return cs;
	}

	public CellStyle bodyStyle() {
		CellStyle cs1 = book.createCellStyle();
		DataFormat format = book.createDataFormat();

		cs1.setAlignment(CellStyle.ALIGN_CENTER);
		cs1.setVerticalAlignment(CellStyle.VERTICAL_CENTER);
		cs1.setBorderLeft(CellStyle.BORDER_THIN);
		cs1.setBorderTop(CellStyle.BORDER_THIN);
		cs1.setBorderRight(CellStyle.BORDER_THIN);
		cs1.setBorderBottom(CellStyle.BORDER_THIN);
		cs1.setDataFormat(format.getFormat("0.0000"));

		return cs1;
	}

}
