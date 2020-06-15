/*  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *  All we ask is that proper credit is given for our work, which includes
 *  - but is not limited to - adding the above copyright notice to the beginning
 *  of your source code files, and to any copyright notice that you may distribute
 *  with programs based on this work.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */
package org.openscience.cdk.similarity;


import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;

import java.util.BitSet;

/**
 *  Calculates the Tanimoto coefficient for a given pair of two 
 *  fingerprint bitsets or real valued feature vectors.
 *
 *  The Tanimoto coefficient is one way to 
 *  quantitatively measure the "distance" or similarity of 
 *  two chemical structures. 
 *
 *  <p>You can use the FingerPrinter class to retrieve two fingerprint bitsets.
 *  We assume that you have two structures stored in cdk.Molecule objects.
 *  A tanimoto coefficient can then be calculated like:
 *  <pre>
 *   BitSet fingerprint1 = Fingerprinter.getFingerprint(molecule1);
 *   BitSet fingerprint2 = Fingerprinter.getFingerprint(molecule2);
 *   float tanimoto_coefficient = Tanimoto.calculate(fingerprint1, fingerprint2);
 *  </pre>
 *
 *  <p>The FingerPrinter assumes that hydrogens are explicitely given, if this 
 *  is desired! 
 *  <p>Note that the continuous Tanimoto coefficient does not lead to a metric space
 *
 *@author         steinbeck
 * @cdk.githash
 *@cdk.created    2005-10-19
 *@cdk.keyword    jaccard
 *@cdk.keyword    similarity, tanimoto
 * @cdk.module fingerprint
 */
@TestClass("org.openscience.cdk.similarity.TanimotoTest")
public class Tanimoto 
{

    /**
     * Evaluates Tanimoto coefficient for two bit sets.
     *
     * @param bitset1 A bitset (such as a fingerprint) for the first molecule
     * @param bitset2 A bitset (such as a fingerprint) for the second molecule
     * @return The Tanimoto coefficient
     * @throws org.openscience.cdk.exception.CDKException  if bitsets are not of the same length
     */
    @TestMethod("testTanimoto1,testTanimoto2")
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
     * Evaluates the continuous Tanimoto coefficient for two real valued vectors.
     *
     * @param features1 The first feature vector
     * @param features2 The second feature vector
     * @return The continuous Tanimoto coefficient
     * @throws org.openscience.cdk.exception.CDKException  if the features are not of the same length
     */
    @TestMethod("testTanimoto3")
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
            c += features1[i] * features2[i];
           // c += Math.abs(features1[i]) * Math.abs(features2[i]);
            a += features1[i]*features1[i];
            b += features2[i]*features2[i];
    	    //c += Math.abs(features1[i] - features2[i]);
    	    //a += Math.abs(features1[i]) + Math.abs(features2[i]);
        }
       
       return c/(a+b-c);
    }
    
    public static double calculate(double features1, double features2) throws CDKException {
        
        // c es la suma del producto de las propiedades de ambos vectores en cada posicion
        // a es la suma del producto de las propiedades del vector 1
        // b es la suma del producto de las propiedades del vector 2
        
        double c = 0.0;
        double a = 0.0;
        double b = 0.0;
        
        c += features1 * features2;
        a += features1 * features1;
        b += features2 * features2 ;
        
        //c += Math.abs(features1 - features2);
  	    //a += Math.abs(features1) + Math.abs(features2);
        
  	    return c/(a+b-c);
  	    //return 1-(c/a);
    }
}
