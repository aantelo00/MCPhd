package org.openscience.cdk.similarity;

public class Squared_Euclidean {
	
	public double calculate(double[] features1, double[] features2) /*throws CDKException*/ {

        /*if (features1.length != features2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }*/

        int n = features1.length;
        
        // c es la suma del producto de las propiedades de ambos vectores en cada posicion
        // a es la suma del producto de las propiedades del vector 1
        // b es la suma del producto de las propiedades del vector 2
        
        double c = 0.0;
        double a = 0.0;
        double b = 0.0;

        for (int i = 0; i < n; i++) {
            c += features1[i] * features2[i];
            a += features1[i]*features1[i];
            b += features2[i]*features2[i];
        }
        
        return ((a+b)-(2*c))/n;
    }

}
