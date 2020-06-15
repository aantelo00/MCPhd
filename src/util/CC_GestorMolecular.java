package util;

import java.util.ArrayList;
import java.util.StringTokenizer;

import domain.CE_Opcion;



public class CC_GestorMolecular {

	private static ArrayList<String> MolSources;
	private static ArrayList<CE_Opcion> Opciones;
	private boolean ismol;

	public CC_GestorMolecular() {
		MolSources = new ArrayList<String>();
		Opciones = new ArrayList<CE_Opcion>();
		ismol = false;
	}

	public String GetMol(int pos)
	{
		return MolSources.get(pos);
	}

	public int GetMolCount()
	{
		return MolSources.size();
	}

	public CE_Opcion GetOpcion(int pos)
	{
		return Opciones.get(pos);
	}
	
	public CE_Opcion GetOpcion(String molName)
	{
		int pos = -1;
		for (int i=0; i<MolSources.size(); i++)
			if (MolSources.get(i).endsWith(molName))
			{
				pos = i;
				break;
			}
		return Opciones.get(pos);
	}
	
	public int GetMolIndex(String molName)
	{
		int pos = -1;
		for (int i=0; i<MolSources.size(); i++)
			if (MolSources.get(i).endsWith(molName))
			{
				pos = i;
				break;
			}
		return pos;
	}

	public void CargarMolecula(String path) {
		if (!MolSources.contains(new String(path)))
		{
			MolSources.add(path);
			CE_Opcion newOpcion = new CE_Opcion();
			ismol  = !path.toLowerCase().endsWith("pdb");
			newOpcion.setMol(ismol);
			Opciones.add(newOpcion);
		}
	}

	public float GetMM(int posicion) {
		return 0;
	}

	public float[] CalculaPC(int posicion) {
		return null;
	}

	public void EliminarMolecula(int posicion) {
		MolSources.remove(posicion);
	}

	public ArrayList<String> getMols() {
		return MolSources;
	}

	public String ParseOpcion(CE_Opcion opcion) {
		/*String StrScript = "select "+opcion.getSelection()
		+"; cpk "+opcion.getAtomCPK()
		+"; wireframe "+opcion.getWireframe()
		+"; label "+opcion.getLabel()
		+"; vector "+opcion.getVector()
		+"; set measure "+opcion.getMeasures()
		+"; zoom "+opcion.getZoom()
		+"; anim "+opcion.getAnim()
		+"; vibration "+opcion.getVibration()
		+"; set selectionhalos "+opcion.isSelectionHalos()
		+"; set showhydrogens "+opcion.isHydShow()
		+"; set showMeasurements "+opcion.isMeasShow()
		+"; set labeloffset "+opcion.getLabelOffset()
		+"; hide "+opcion.getHide()
		+"; set PerspectiveDepth "+opcion.isPerspective()
		+"; set showBoundBox "+opcion.isBoundBox()
		+"; set showAxes "+opcion.isAxes()
		+"; "+opcion.getStructure();*/
		return new String();//StrScript;
	}

	public void UpdatePos(int selectedMolPos, CE_Opcion op) {
		ArrayList<CE_Opcion> temp = new ArrayList<CE_Opcion>();
		for (int i=0; i<Opciones.size(); i++)
			if (i==selectedMolPos)
				temp.add(op);
			else temp.add(Opciones.get(i));
		Opciones = temp;
	}

	public void EliminarOpcion(int i) {
		Opciones.remove(i);
	}

	public String getNoDefaultOptions(int selectionIndex) {
		String retorno = "";
		CE_Opcion c = new CE_Opcion();
		CE_Opcion c1 = Opciones.get(selectionIndex);
		if (c1.getSelection().compareToIgnoreCase(c.getSelection()) != 0)	retorno += "select "+c1.getSelection()+"\n";
		if (c1.getAnim().compareToIgnoreCase(c.getAnim()) != 0)				retorno += "anim "+c1.getAnim()+"\n";
		if (c1.getAtomCPK().compareToIgnoreCase(c.getAtomCPK()) != 0)		retorno += "cpk "+c1.getAtomCPK()+"\n";
		if (c1.getHide().compareToIgnoreCase(c.getHide()) != 0)				retorno += "hide "+c1.getHide()+"\n";
		if (c1.getWireframe().compareToIgnoreCase(c.getWireframe()) != 0)	retorno += "wireframe "+c1.getWireframe()+"\n";
		if (c1.getLabel().compareToIgnoreCase(c.getLabel()) != 0)			retorno += "label "+c1.getLabel()+"\n";
		if (c1.getVector().compareToIgnoreCase(c.getVector()) != 0)			retorno += "vector "+c1.getVector()+"\n";
		if (c1.getMeasures().compareToIgnoreCase(c.getMeasures()) != 0)		retorno += "set measure "+c1.getMeasures()+"\n";
		if (c1.getZoom().compareToIgnoreCase(c.getZoom()) != 0)				retorno += "zoom "+c1.getZoom()+"\n";
		if (c1.getVibration().compareToIgnoreCase(c.getVibration()) != 0)	retorno += "vibration "+c1.getVibration()+"\n";
		if (c1.isSelectionHalos() != c.isSelectionHalos())					retorno += "set selectionhalos "+c1.isSelectionHalos()+"\n";
		if (c1.isHydShow() != c.isHydShow())								retorno += "set showhydrogens "+c1.isHydShow()+"\n";
		if (c1.isMeasShow() != c.isMeasShow())								retorno += "set showMeasurements "+c1.isMeasShow()+"\n";
		if (c1.getLabelOffset().compareToIgnoreCase(c.getLabelOffset()) != 0)	retorno += "set labeloffset "+c1.getLabelOffset()+"\n";
		if (c1.isPerspective() != c.isPerspective())						retorno += "set PerspectiveDepth "+c1.isPerspective()+"\n";
		if (c1.isBoundBox() != c.isBoundBox())								retorno += "set showBoundBox "+c1.isBoundBox()+"\n";
		if (c1.isAxes() != c.isAxes())										retorno += "set showAxes "+c1.isAxes()+"\n";
		if (c1.getStereo().compareToIgnoreCase(c.getStereo()) != 0)			retorno += "stereo "+c1.getLabelOffset()+"\n";
		if (c1.isRotating()) 												retorno += "spin on\n";
		return retorno;
	}

	public void ActualizaOpcion(CE_Opcion op, String script) {
		StringTokenizer t1 = new StringTokenizer(script);
		while (t1.hasMoreTokens())
		{
			String token = t1.nextToken("\n");
			if (token.contains("all") || token.contains("none") || token.contains("hydrogen") || token.contains("carbon") || token.contains("nitrogen") || token.contains("oxygen")
					|| token.contains("phosphorus") || token.contains("sulphur") || token.contains("fluorine") || token.contains("chlorine") || token.contains("iodine") || token.contains("bromine")){
				op.setSelection(token.substring(token.lastIndexOf(' ') + 1,token.length()));
			}
			if (token.contains("anim")) 									op.setAnim(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("cpk")) 										op.setAtomCPK(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("hide")) 									op.setHide(token.substring(token.indexOf(' ') + 1));
			if (token.contains("wireframe"))								op.setWireframe(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("label")) 									op.setLabel(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("vector")) 									op.setVector(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("measure")) 									op.setMeasures(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("zoom")) 									op.setZoom(token.substring(token.lastIndexOf(' ')+1));
			if (token.contains("vibration"))								op.setVibration(token.substring(token.lastIndexOf("\n")+1));
			if (token.contains("showhydrogens")) 							op.setHydShow(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("showMeasurements"))							op.setMeasShow(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("labeloffset")) 								op.setLabelOffset(token.substring(token.lastIndexOf("labeloffset")+1));
			if (token.contains("PerspectiveDepth")) 						op.setPerspective(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("showBoundBox")) 							op.setBoundBox(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("showAxes")) 								op.setAxes(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("stereo")) 									op.setStereo(token.substring(token.indexOf(' ')+1));
			if (token.contains("isosurface") || token.contains("dots")) 	op.setSurface(token);
			if (token.contains("on")){
				op.setRotating(true);
			}else{
				op.setRotating(false);
			}
			if (token.contains("selectionhalos"))   						op.setSelectionHalos(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
			if (token.contains("showBoundBox"))   							op.setBoundBox(Boolean.parseBoolean(token.substring(token.lastIndexOf(' ')+1)));
		}
	}

	public int GetOpcionCount() {
		return Opciones.size();
	}

	public void AddOpcion(CE_Opcion opcion) {
		Opciones.add(opcion);
	}

	public void SetMoleculePropertiesCalulated(int selectedMolPos, boolean selected) {
		for (int i=0; i<Opciones.size(); i++)
			if (i == selectedMolPos)
				Opciones.get(i).setCalculos(selected);
			else Opciones.get(i).setCalculos(false);

	}

	public String getDefaultOption() {
		String retorno = "";
		CE_Opcion c1 = new CE_Opcion();
		retorno += "select "+c1.getSelection()+"\n";
		retorno += "anim "+c1.getAnim()+"\n";
		retorno += "cpk "+c1.getAtomCPK()+"\n";
		retorno += "hide "+c1.getHide()+"\n";
		retorno += "wireframe "+c1.getWireframe()+"\n";
		retorno += "label "+c1.getLabel()+"\n";
		retorno += "vector "+c1.getVector()+"\n";
		retorno += "set measure "+c1.getMeasures()+"\n";
		retorno += "zoom "+c1.getZoom()+"\n";
		retorno += "vibration "+c1.getVibration()+"\n";
		retorno += "set selectionhalos "+c1.isSelectionHalos()+"\n";
		retorno += "set showhydrogens "+c1.isHydShow()+"\n";
		retorno += "set showMeasurements "+c1.isMeasShow()+"\n";
		retorno += "set labeloffset "+c1.getLabelOffset()+"\n";
		retorno += "set PerspectiveDepth "+c1.isPerspective()+"\n";
		retorno += "set showBoundBox "+c1.isBoundBox()+"\n";
		retorno += "set showAxes "+c1.isAxes()+"\n";
		if (c1.isRotating()) retorno += "spin on\n";
		return retorno;
	}
}
