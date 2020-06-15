package org.openscience.cdk.similarity;

public class Ruzicka {
	
	public double calculate(double[] features1, double[] features2) /*throws CDKException*/ {

        /*if (features1.length != features2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }*/

		int n = features1.length;
        double miab = 0.0;
        double maab = 0.0;

        for (int i = 0; i < n; i++) {
            maab += Math.max(features1[i], features2[i]);
            miab += Math.min(features1[i], features2[i]);

        }
        return miab / maab;
    }

}
