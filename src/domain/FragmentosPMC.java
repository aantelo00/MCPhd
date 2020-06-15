package domain;

public class FragmentosPMC {
	private FragmentoMolecular fragMol1;
	private FragmentoMolecular fragMol2;
	private double indexSimilary;
	
	public FragmentosPMC(FragmentoMolecular fragMol1, FragmentoMolecular fragMol2) {
		//super();
		this.fragMol1 = fragMol1;
		this.fragMol2 = fragMol2;
		indexSimilary = 0.0;
	}

	/**
	 * @return the indexSimilary
	 */
	public double getIndexSimilary() {
		return indexSimilary;
	}

	/**
	 * @param indexSimilary the indexSimilary to set
	 */
	public void setIndexSimilary(double indexSimilary) {
		this.indexSimilary = indexSimilary;
	}

	/**
	 * @return the fragMol1
	 */
	public FragmentoMolecular getFragMol1() {
		return fragMol1;
	}

	/**
	 * @param fragMol1 the fragMol1 to set
	 */
	public void setFragMol1(FragmentoMolecular fragMol1) {
		this.fragMol1 = fragMol1;
	}

	/**
	 * @return the fragMol2
	 */
	public FragmentoMolecular getFragMol2() {
		return fragMol2;
	}

	/**
	 * @param fragMol2 the fragMol2 to set
	 */
	public void setFragMol2(FragmentoMolecular fragMol2) {
		this.fragMol2 = fragMol2;
	}
	
	
	
	

}
