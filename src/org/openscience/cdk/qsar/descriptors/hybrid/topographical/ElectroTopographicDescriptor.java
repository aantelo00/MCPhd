package org.openscience.cdk.qsar.descriptors.hybrid.topographical;

import java.util.Iterator;
import java.util.List;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.DescriptorSpecification;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IAtomicDescriptor;
import org.openscience.cdk.qsar.IHybridDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.DistanceToAtomDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.PeriodicTable;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class ElectroTopographicDescriptor implements IHybridDescriptor {

	private IAtomContainer container;
	
	@Override
	public DescriptorValue calculate(IAtomContainer atomcontainer) {

		try {
			container = (IAtomContainer) atomcontainer.clone();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);            
			CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(container.getBuilder());
			hAdder.addImplicitHydrogens(container);
			AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
			boolean aromatic = CDKHueckelAromaticityDetector.detectAromaticity(container);
			container.setFlag(CDKConstants.ISAROMATIC, aromatic);
		} catch (CloneNotSupportedException e) {
			return getDummyDescriptorValue(new CDKException("Error during clone"));
		} catch (CDKException e) {
			return getDummyDescriptorValue(new CDKException("Error during atom typing" + e.getMessage()));
		}

		if(!isEIntrinsicValue()){
			
			for(IAtom atom : this.container.atoms()){
				String symbol = atom.getSymbol();
				if(!symbol.equals("H")){
					int periodo = PeriodicTable.getPeriod(symbol);
					int electronesValencia = atom.getValency();
					int delta = atom.getFormalNeighbourCount() - numeroHidrogenos(atom);
					double estadoElecTopologico = (Math.pow(2/periodo,2)* electronesValencia + 1)/delta;
					atom.setProperty(IHybridDescriptor.Descriptor.ELECTROINTRINSICVALUE, estadoElecTopologico);
				}
			}
		}
		
		double total = 0.0;
		try {

			total = updateRefractedTopographicValues();
			this.container.setProperty(IHybridDescriptor.Descriptor.TOTALELECTROTOPOGRAPHIC, total);

		} catch (CDKException e) {
			return getDummyDescriptorValue(new CDKException("Error during 3D coordinates" + e.getMessage()));
		}

		return new DescriptorValue(
				getSpecification(), getParameterNames(), getParameters(), 
				new DoubleResult(total),getDescriptorNames());

	}

	private double updateRefractedTopographicValues() throws CDKException {

		IAtomicDescriptor distance = new DistanceToAtomDescriptor();
		double efecto = 0.0, total = 0.0;
		//Object[] arrAux = new Object[1];

		int nAtomos = container.getAtomCount();
		for (int i = 0; i < nAtomos; i++) {
			efecto = 0.0;
			IAtom atomoOrigen = container.getAtom(i);
			if(!atomoOrigen.getSymbol().equals("H")){
				double valorEstateIntAtomoOrigen = Double.valueOf(atomoOrigen.getProperty(IHybridDescriptor.Descriptor.ELECTROINTRINSICVALUE).toString());
				for(int j = 0; j < nAtomos; j++){
					IAtom atomoDestino = container.getAtom(j);
					if(i != j && !atomoDestino.getSymbol().equals("H")){

						//arrAux[0] = j; 
						Object[] objs = {j};
						distance.setParameters(objs);
						//(IntegerArrayResult)descriptor.calculate(containersList.get(0)).getValue();
						DescriptorValue value =  distance.calculate(atomoOrigen, container);
						double dist = ((DoubleResult) value.getValue()).doubleValue();
						if(dist != 0){
							double valorEstateIntAtomoDestino = Double.valueOf(atomoDestino.getProperty(IHybridDescriptor.Descriptor.ELECTROINTRINSICVALUE).toString());
							efecto += (valorEstateIntAtomoOrigen - valorEstateIntAtomoDestino)/Math.pow(dist, 2);

						}
					}
				}

				atomoOrigen.setProperty(IHybridDescriptor.Descriptor.ELECTROTOPOGRAPHIC, valorEstateIntAtomoOrigen + efecto);
				total += valorEstateIntAtomoOrigen + efecto;
			}
		}

		return total;

	}

	private boolean isEIntrinsicValue(){

		boolean aux = false;
		try{
			//Primer atomo, si esta bien formado no debe ser H
			this.container.getFirstAtom().getProperty(IHybridDescriptor.Descriptor.ELECTROINTRINSICVALUE).toString();
			aux = true;

		}catch (Exception e) {

		}

		return aux;
	}

	private int numeroHidrogenos(IAtom atom) {

		int cHidrogeno = 0;
		List<IAtom> atomosConectados = this.container.getConnectedAtomsList(atom);
		for (Iterator<IAtom> atomo = atomosConectados.iterator(); atomo.hasNext();) {
			IAtom iAtom = (IAtom) atomo.next();
			if (iAtom.getSymbol().equals("H"))
				cHidrogeno++;
		}

		return cHidrogeno;
	}

	private DescriptorValue getDummyDescriptorValue(Exception e) {
		DoubleArrayResult results = new DoubleArrayResult();
		results.add(Double.NaN);
		results.add(Double.NaN);
		results.add(Double.NaN);
		return new DescriptorValue(getSpecification(), getParameterNames(),
				getParameters(), results, getDescriptorNames(), e);
	}
	
	@Override
	public IAtomContainer getAtomContainer() {

		return container;
	}

	@Override
	public String[] getDescriptorNames() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String[] getParameterNames() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object getParameterType(String name) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Object[] getParameters() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public DescriptorSpecification getSpecification() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setParameters(Object[] params) throws CDKException {
		// TODO Auto-generated method stub

	}

}
