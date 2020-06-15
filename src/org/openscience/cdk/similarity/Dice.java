package org.openscience.cdk.similarity;


import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;

import java.util.BitSet;


@TestClass("org.openscience.cdk.similarity.TanimotoTest")
public class Dice
{
	

    /**
     * Evaluates Coseno coefficient for two bit sets.
     *
     * @param bitset1 A bitset (such as a fingerprint) for the first molecule
     * @param bitset2 A bitset (such as a fingerprint) for the second molecule
     * @return The Coseno coefficient
     * @throws org.openscience.cdk.exception.CDKException  if bitsets are not of the same length
     */    
    @TestMethod("testCoseno1,testCoseno2")
    public static float calculate(BitSet bitset1, BitSet bitset2) throws CDKException
    {
        float _bitset1_cardinality = bitset1.cardinality();
        float _bitset2_cardinality = bitset2.cardinality();
        if (bitset1.size() != bitset2.size()) {
            throw new CDKException("Bisets must have the same bit length");
        }
        BitSet one_and_two = (BitSet)bitset1.clone();
        one_and_two.and(bitset2);
        float _common_bit_count = one_and_two.cardinality();
        return _common_bit_count/(_bitset1_cardinality + _bitset2_cardinality - _common_bit_count);
    }
    
    /**
     * Evaluates the continuous Coseno coefficient for two real valued vectors.
     *
     * @param features1 The first feature vector
     * @param features2 The second feature vector
     * @return The continuous Coseno coefficient
     * @throws org.openscience.cdk.exception.CDKException  if the features are not of the same length
     */
    @TestMethod("testCoseno3")
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
            c += Math.abs(features1[i] * features2[i]);
            a += features1[i]*features1[i];
            b += features2[i]*features2[i];
        }
        return c / ((a+b)*0.5);
    }
}
