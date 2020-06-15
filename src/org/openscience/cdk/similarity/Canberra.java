package org.openscience.cdk.similarity;

import org.openscience.cdk.exception.CDKException;

public class Canberra {
	
	public static double calculate(double[] features1, double[] features2) throws CDKException {
		 
    	if (features1.length != features2.length) {
            throw new CDKException("Features vectors must be of the same length");
        }

        int n = features1.length;
        
        // c es la suma del producto de las propiedades de ambos vectores en cada posicion
        // a es la suma del producto de las propiedades del vector 1
        // b es la suma del producto de las propiedades del vector 2
        
        double c = 0.0;
        double a = 0.0;
        double b = 0.0;

      for (int i = 0; i < n; i++) {
    	    a += Math.abs(features1[i] - features2[i]);
    	    b += Math.abs(features1[i]) + Math.abs(features2[i]);
    	    c+= a/b;
        }
       
       return c;
    }
    
    public static double calculate(double features1, double features2) throws CDKException {
        
        // c es la suma del producto de las propiedades de ambos vectores en cada posicion
        // a es la suma del producto de las propiedades del vector 1
        // b es la suma del producto de las propiedades del vector 2
        
        double a = 0.0;
        double b = 0.0;
        
        a += Math.abs(features1 - features2);
  	    b += Math.abs(features1) + Math.abs(features2);
        
  	    return a/b;
    }

}
