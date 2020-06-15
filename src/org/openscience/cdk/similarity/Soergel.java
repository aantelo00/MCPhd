package org.openscience.cdk.similarity;

public class Soergel {

    public double calculate(double[] features1, double[] features2) /*throws CDKException*/ {

        /*if (features1.length != features2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }*/

    	int n = features1.length;
        double maab = 0.0;
        double r = 0.0;
        for (int i = 0; i < n; i++) {
            r += Math.abs(features1[i] - features2[i]);
            maab += Math.max(features1[i], features2[i]);
        }
        return r / maab;
    }

}
