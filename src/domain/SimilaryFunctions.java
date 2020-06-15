package domain;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.similarity.Canberra;
//import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.similarity.Cosine;
import org.openscience.cdk.similarity.Dice;
/*import org.openscience.cdk.similarity.Cosine;
import org.openscience.cdk.similarity.Squared_Euclidean;
import org.openscience.cdk.similarity.Russell_Rad;
import org.openscience.cdk.similarity.Forbes;
import org.openscience.cdk.similarity.Simpson;
import org.openscience.cdk.similarity.Pearson;
import org.openscience.cdk.similarity.Dennis;
import org.openscience.cdk.similarity.Sokal_Sneath;
import org.openscience.cdk.similarity.Kulczynski;
import org.openscience.cdk.similarity.Ochiai;
import org.openscience.cdk.similarity.Czekanowski;*/
//import org.openscience.cdk.similarity.Ruzicka;
/*import org.openscience.cdk.similarity.Soergel;
import org.openscience.cdk.similarity.S�rensen;*/
import org.openscience.cdk.similarity.Tanimoto;
/*import org.openscience.cdk.similarity.VicisWave;*/





import domain.TypeSimilaryFunction;

public class SimilaryFunctions {

	public SimilaryFunctions() {

	}
	
	@SuppressWarnings("static-access")
	public double CalculeSimilaryFunctionAtomo(TypeSimilaryFunction type, double[] atom1, double[] atom2) throws CDKException {
		double value = 0;

		switch (type) {
			case TANIMOTO: {
				value = new Tanimoto().calculate(atom1, atom2);
				break;
			}
		default:
			break;
		}
		return value;
	}

	@SuppressWarnings("static-access")
	public double CalculeSimilaryFunction(TypeSimilaryFunction type, double[] baseFragment, double[] Fragment) throws CDKException {
		double value = 0;
		/*switch (descriptor) {
		  case 'E': {
			  descriptorBaseFrament = new double[1]; 
			  descriptorFragment = new double[1];
			  descriptorBaseFrament[0] = baseFragment[0];
			  descriptorFragment[0] = Fragment[0];
		  }
		  case 'R': {
			  descriptorBaseFrament = new double[1]; 
			  descriptorFragment = new double[1];
			  descriptorBaseFrament[1] = baseFragment[1];
			  descriptorFragment[1]= Fragment[1];
		  }
		  case 'L': {
			  descriptorBaseFrament = new double[1]; 
			  descriptorFragment = new double[1];
			  descriptorBaseFrament[2] = baseFragment[2];
			  descriptorFragment[2] = Fragment[2];
		  }
		  case 'A': {
			  descriptorBaseFrament = new double[3]; 
			  descriptorFragment = new double[3];
			  descriptorBaseFrament = baseFragment;
			  descriptorFragment = Fragment;
		  }
		}*/

		switch (type) {
			case TANIMOTO: {
				value = new Tanimoto().calculate(baseFragment, Fragment);
				break;
			}
			/*case COSINE: {
				value = new Cosine().calculate(baseFragment, Fragment);
				break;
			}
			case DICE: {
				value = new Dice().calculate(baseFragment, Fragment);
				break;
			}*/
			/*case SQUARED_EUCLIDEAN: {
				value = new Squared_Euclidean().calculate(baseFragment, Fragment);
				break;
			}
			case RUSSELL_RAD: {
				value = new Russell_Rad().calculate(baseFragment, Fragment);
				break;
			}
			case FORBES: {
				value = new Forbes().calculate(baseFragment, Fragment);
				break;
			}*/
			/*case SIMPSON: {
				value = new Simpson().calculate(baseFragment, Fragment);
				break;
			}*/
			/*case PEARSON: {
				value = new Pearson().calculate(baseFragment, Fragment);
				break;
			}
			case DENNIS: {
				value = new Dennis().calculate(baseFragment, Fragment);
				break;
			}
			case SOKAL_SNEATH: {
				value = new Sokal_Sneath().calculate(baseFragment, Fragment);
				break;
			}
			case KULCZYNSKI: {
				value = new Kulczynski().calculate(baseFragment, Fragment);
				break;
			}
			case OCHIAI: {
				value = new Ochiai().calculate(baseFragment, Fragment);
				break;
			}
			case SORENSEN: {
				value = new S�rensen().calculate(baseFragment, Fragment);
				break;
			}
			case SOERGEL: {
				value = new Soergel().calculate(baseFragment, Fragment);
				break;
			}
			case CZEKANOWSKI: {
				value = new Czekanowski().calculate(baseFragment, Fragment);
				break;
			}*/
			/*case RUZICKA: {
				value = new Ruzicka().calculate(baseFragment, Fragment);
				break;
			}*/
			/*case VICIS_WAVE: {
				value = new VicisWave().calculate(baseFragment, Fragment);
				break;
			}*/
			
		}
		return value;
	}
	
	@SuppressWarnings("static-access")
	public double CalculeSimilaryFunction(TypeSimilaryFunction type, double baseFragment, double Fragment) throws CDKException {
		double value = 0;

		switch (type) {
			case TANIMOTO: {
				value = new Tanimoto().calculate(baseFragment, Fragment);
				break;
			}
			/*case CANBERRA: {
				value = 1- new Canberra().calculate(baseFragment, Fragment);
				//System.out.println(baseFragment + " - " + Fragment + " - " + value);
				break;
			}*/
		default:
			break;
	
		}
		return value;
	}
}
