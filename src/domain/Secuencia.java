package domain;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

import domain.TipoCentroDescriptor;

public class Secuencia {

	protected StringBuffer stb;
	public String getCodigo(IAtomContainer ac, TipoCentroDescriptor tipo){
		String secuencia = "";
		if(tipo == TipoCentroDescriptor.ANILLO || tipo == TipoCentroDescriptor.CLUSTER3 || tipo == TipoCentroDescriptor.CLUSTER4)
		{
			for (int i = 0; i < ac.getAtomCount(); i++) {
				if(i < ac.getAtomCount()-1)
					secuencia += ac.getAtom(i).getSymbol() + ac.getAtom(i).getID() + ",";
				else secuencia += ac.getAtom(i).getSymbol() + ac.getAtom(i).getID();
					
			}
		} else secuencia += ac.getAtom(0).getSymbol() + ac.getAtom(0).getID();
		
		return secuencia;
		
	}
	
	protected String codificarAtomo(IBond bond, IAtom at2){

		if(bond.getFlag(CDKConstants.ISAROMATIC)) {
		} else{
			switch (bond.getOrder()) {
			case SINGLE:
				break;
			case DOUBLE:
				break;
			case TRIPLE:
				break;
			default: 
				break;
			}
		}

		//String codigos = orden + at2.getSymbol() + at2.getID();
		String codigos = "," + at2.getSymbol() + at2.getID() ;
		return codigos;
	}
	
}
