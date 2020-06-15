package org.openscience.cdk.similarity;

public class VicisWave {
	
	public double calculate(double[] numeralSet1, double[] numeralSet2) /*throws CDKException*/ {

        /*if (features1.length != features2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }*/

	    int quantity = numeralSet1.length;
        double module = 0.0;
        double min = 0.0;
        double value = 0.0;

        for (int i = 0; i < quantity; i++) {
        	module = Math.pow(Math.abs(numeralSet1[i] - numeralSet2[i]),2);
            min = Math.pow(Math.min(numeralSet1[i], numeralSet2[i]),2);
            value += module / min;
        }
        return Math.abs(value);

	} 
}
	
