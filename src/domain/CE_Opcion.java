package domain;

import java.util.ArrayList;

public class CE_Opcion {

	private String AtomCPK, Wireframe, Label, Vector, Measures, Zoom, Anim, Vibration, Hbonds, SBonds, LabelOffset, Hide, mainV, Ac, protein, hetero;
	private boolean Calculos;
	private boolean Geometria;
	private String selection;
	private boolean SelectionHalos, HydShow, MeasShow, Perspective, Axes, BoundBox, Mol, Rotating;
	private String Stereo;
	private String Surface, Structure, StructureOpacity,stateSelection,stateAtom;
	private ArrayList<String> listaAtomos;
	private ArrayList<String> listaAtomosState;
				
	public CE_Opcion() {
		selection = "all";
		Rotating = false;
		AtomCPK = "20%";
		Wireframe = ".15";
		Label = "off";
		Vector = "off";
		Measures = "nanometers";
		Zoom = "100";
		Anim = "off";
		Vibration = "off";
		SelectionHalos = false;
		Hbonds = "off";
		SBonds = "off";
		HydShow = true;
		MeasShow = true;
		LabelOffset = "4 4";
		Perspective = true;
		Axes = false;
		BoundBox = false;
		Hide = "none";
		Structure = "backbone off;cartoons off;ribbons off;rockets off;strands off;trace off";
		Stereo = "off";
		Surface = "isosurface delete;select *;dots off";
		protein = "protein";
		hetero = "hetero";
		mainV = "front";
		Ac    = "off";
		Mol = true;
		Calculos = false;
		Geometria = false;
		StructureOpacity = "translucent";
		stateSelection = "";
		stateAtom = "20%";
		listaAtomos = new ArrayList<String>();
		listaAtomosState = new ArrayList<String>();
	}
	
	public String getAc() {
		return Ac;
	}

	public String getAnim() {
		return Anim;
	}

	public String getAtomCPK() {
		return AtomCPK;
	}


	public String getHbonds() {
		return Hbonds;
	}


	public String getHetero() {
		return hetero;
	}


	public String getHide() {
		return Hide;
	}


	public String getLabel() {
		return Label;
	}


	public String getLabelOffset() {
		return LabelOffset;
	}
	
	public ArrayList<String> getListaAtomos() {
		return listaAtomos;
	}
	
	public ArrayList<String> getListaAtomosState() {
		return listaAtomosState;
	}

	public String getMainV() {
		return mainV;
	}


	public String getMeasures() {
		return Measures;
	}


	public String getProtein() {
		return protein;
	}


	public String getSBonds() {
		return SBonds;
	}


	public String getSelection() {
		return selection;
	}


	public String getStereo() {
		return Stereo;
	}


	public String getStructure() {
		return Structure;
	}
	
	public String getStructureOpacity() {
		return StructureOpacity;
	}


	public String getSurface() {
		return Surface;
	}

	public String getStateSelection() {
		return stateSelection;
	}
	
	public String getStateAtom() {
		return stateAtom;
	}
	
	public String getVector() {
		return Vector;
	}


	public String getVibration() {
		return Vibration;
	}


	public String getWireframe() {
		return Wireframe;
	}


	public String getZoom() {
		return Zoom;
	}


	public boolean isAxes() {
		return Axes;
	}


	public boolean isBoundBox() {
		return BoundBox;
	}
	
	public boolean isCalculo() {
		return Calculos;
	}


	public boolean isHydShow() {
		return HydShow;
	}


	public boolean isMeasShow() {
		return MeasShow;
	}


	public boolean isMol() {
		return Mol;
	}


	public boolean isOptimizada() {
	    return Geometria;
	}


	public boolean isPerspective() {
		return Perspective;
	}


	public boolean isRotating() {
		return Rotating;
	}


	public boolean isSelectionHalos() {
		return SelectionHalos;
	}


	public void RestoreDefaultStructureAndSurface()
	{
		Structure = "backbone off;cartoons off;ribbons off;rockets off;strands off;trace off";
		Surface   = "isosurface delete;select *;dots off";
	}


	public void setAc(String ac) {
		Ac = ac;
	}


	public void setAnim(String anim) {
		Anim = anim;
	}


	public void setAtomCPK(String atomCPK) {
		AtomCPK = atomCPK;
	}


	public void setAxes(boolean axes) {
		Axes = axes;
	}

	public void setBoundBox(boolean box) {
		BoundBox = box;
		
	}

	public void setCalculos(boolean b) {
		Calculos = b;
		
	}

	public void setGeometria(boolean b) {
		Geometria = b;
		
	}

	public void setHbonds(String hbonds) {
		Hbonds = hbonds;
	}

	public void setHetero(String hetero) {
		this.hetero = hetero;
	}

	public void setHide(String hide) {
		Hide = hide;
	}

	public void setHydShow(boolean hydshow) {
		HydShow = hydshow;
	}

	public void setLabel(String label) {
		Label = label;
	}

	public void setLabelOffset(String string) {
		LabelOffset = string;
	}
	
	public void setListaAtomos(ArrayList<String> listaAtomos) {
		this.listaAtomos = listaAtomos;
	}
	
	public void setListaAtomosState(ArrayList<String> listaAtomosState) {
		this.listaAtomosState = listaAtomosState;
	}
	
	public void setMainV(String mainV) {
		this.mainV = mainV;
	}

	public void setMeasShow(boolean measShow) {
		MeasShow = measShow;
	}

	public void setMeasures(String measures) {
		Measures = measures;
	}
	
	public void setMol(boolean mol) {
		Mol = mol;
	}

	public void setPerspective(boolean persp) {
		Perspective = persp;
	}

	public void setProtein(String protein) {
		this.protein = protein;
	}

	public void setRotating(boolean rotating) {
		Rotating = rotating;
	}
	
	public void setSBonds(String bonds) {
		SBonds = bonds;
	}

	public void setSelection(String select) {
		selection = select;
	}
	
	public void setStateSelection(String stateSelection) {
		this.stateSelection = stateSelection;
	}
	
	public void setStateAtom(String stateAtom) {
		this.stateAtom = stateAtom;
	}	
	
	public void setSelectionHalos(boolean selectionHalos) {
		SelectionHalos = selectionHalos;
	}

	public void setStereo(String string) {
		Stereo = string;
	}
	
	public void setStructure(String string) {
		Structure = string;		
	}


	public void setStructureOpacity(String structureOpacity) {
		StructureOpacity = structureOpacity;
	}


	public void setSurface(String string) {
		Surface = string;		
	}


	public void setVector(String vector) {
		Vector = vector;
	}


	public void setVibration(String vibration) {
		Vibration = vibration;
	}

	public void setWireframe(String wireframe) {
		Wireframe = wireframe;
	}

	public void setZoom(String zoom) {
		Zoom = zoom;
	}
	
	public String getNoDefaultOptions() {
		String retorno = "";
		CE_Opcion c = new CE_Opcion();
		if (getSelection().compareToIgnoreCase(c.getSelection()) != 0)		retorno += "select "+getSelection()+"\n";
		if (getAnim().compareToIgnoreCase(c.getAnim()) != 0)				retorno += "anim "+getAnim()+"\n";
		if (getAtomCPK().compareToIgnoreCase(c.getAtomCPK()) != 0)			retorno += "cpk "+getAtomCPK()+"\n";
		if (getHide().compareToIgnoreCase(c.getHide()) != 0)				retorno += "hide "+getHide()+"\n";
		if (getWireframe().compareToIgnoreCase(c.getWireframe()) != 0)		retorno += "wireframe "+getWireframe()+"\n";
		if (getLabel().compareToIgnoreCase(c.getLabel()) != 0)				retorno += "label "+getLabel()+"\n";
		if (getVector().compareToIgnoreCase(c.getVector()) != 0)			retorno += "vector "+getVector()+"\n";
		if (getMeasures().compareToIgnoreCase(c.getMeasures()) != 0)		retorno += "set measure "+getMeasures()+"\n";
		if (getZoom().compareToIgnoreCase(c.getZoom()) != 0)				retorno += "zoom "+getZoom()+"\n";
		if (getVibration().compareToIgnoreCase(c.getVibration()) != 0)		retorno += "vibration "+getVibration()+"\n";
		if (isSelectionHalos() != c.isSelectionHalos())						retorno += "set selectionhalos "+isSelectionHalos()+"\n";
		if (isHydShow() != c.isHydShow())									retorno += "set showhydrogens "+isHydShow()+"\n";
		if (isMeasShow() != c.isMeasShow())									retorno += "set showMeasurements "+isMeasShow()+"\n";
		if (getLabelOffset().compareToIgnoreCase(c.getLabelOffset()) != 0)	retorno += "set labeloffset "+getLabelOffset()+"\n";
		if (isPerspective() != c.isPerspective())							retorno += "set PerspectiveDepth "+isPerspective()+"\n";
		if (isBoundBox() != c.isBoundBox())									retorno += "set showBoundBox "+isBoundBox()+"\n";
		if (isAxes() != c.isAxes())											retorno += "set showAxes "+isAxes()+"\n";
		if (getStereo().compareToIgnoreCase(c.getStereo()) != 0)			retorno += "stereo "+getLabelOffset()+"\n";
		if (isRotating()){
			retorno += "spin on\n";
		}else{
			retorno += "spin off\n";
		}
		return retorno;
	}				
}
