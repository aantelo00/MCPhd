/* $Revision$ $Author$ $Date$ 
 *
 * Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 * (or see http://www.gnu.org/copyleft/lesser.html)
 */
package org.openscience.cdk.isomorphism.matchers.smarts;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IRingSet;

/**
 * This query atom matches any atom with a certain number of SSSR. 
 *
 * @cdk.module  smarts
 * @cdk.githash
 * @cdk.keyword SMARTS 
 */
public class RingMembershipAtom extends SMARTSAtom {
	private static final long serialVersionUID = -7963168231557641862L;

    /** Number of rings to which this atom belongs, if < 0 check any ring membership. */
	private int numSSSR;

    /**
     * Ring membership query atom. Check if the an atom belongs to <i>num</i> of
     * rings. To specify any ring membership, <i>num</i> should be specified as
     * < 0. Generally in SMARTS it's better negate ring membership with {@code
     * [!R]} however for legacy reasons {@code [R0]} was accepted and checks
     * this atoms belongs to 0 rings.
     *
     * @param num number of rings which this atom belongs to, < 0 any ring.
     */
    public RingMembershipAtom(int num) {
		this.numSSSR = num;
	}

	/* (non-Javadoc)
	 * @see org.openscience.cdk.isomorphism.matchers.smarts.SMARTSAtom#matches(org.openscience.cdk.interfaces.IAtom)
	 */
	public boolean matches(IAtom atom) {
		if (atom.getFlag(CDKConstants.ISINRING)) {
			IRingSet ringSet = (IRingSet)atom.getProperty(CDKConstants.SMALLEST_RINGS);
            // < 0 means any ring, as you can see below R0 is valid
			return numSSSR < 0 || ringSet.getAtomContainerCount() == numSSSR;
		} else {
            if (numSSSR == 0) return true;
        }
		return false;
	}
}
