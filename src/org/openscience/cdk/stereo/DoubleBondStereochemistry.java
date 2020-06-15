/* Copyright (C) 2012  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.stereo;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ICDKObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import java.util.Map;

/**
 * Stereochemistry specification for double bonds. See {@link IDoubleBondStereochemistry} for
 * further details.
 *
 * @cdk.module core
 *
 * @see org.openscience.cdk.interfaces.IDoubleBondStereochemistry
 */
@TestClass("org.openscience.cdk.stereo.DoubleBondStereochemistryTest")
public class DoubleBondStereochemistry implements IDoubleBondStereochemistry {

	private Conformation stereo;
	private IBond[] ligandBonds;
	private IBond stereoBond;
    private IChemObjectBuilder builder;

	/**
	 * Creates a new double bond stereo chemistry. The path of length three is defined by
	 * <code>ligandBonds[0]</code>, <code>stereoBonds</code>, and <code>ligandBonds[1]</code>.
	 */
	public DoubleBondStereochemistry(IBond stereoBond, IBond[] ligandBonds, Conformation stereo) {
        if(ligandBonds.length > 2)
            throw new IllegalArgumentException("expected two ligand bonds");
		this.stereoBond = stereoBond;
		this.ligandBonds = ligandBonds;
		this.stereo = stereo;
	}
	
    /**
     * Sets a new {@link IChemObjectBuilder}.
     *
     * @param builder the new {@link IChemObjectBuilder} to be returned
     * @see #getBuilder()
     */
    @TestMethod("testBuilder")
    public void setBuilder(IChemObjectBuilder builder) {
        this.builder = builder;
    }

	/** {@inheritDoc} */
	@Override
    @TestMethod("testBuilder")
	public IChemObjectBuilder getBuilder() {
		return this.builder;
	}

	/** {@inheritDoc} */
	@Override
	@TestMethod("testGetBonds")
	public IBond[] getBonds() {
		IBond[] arrayCopy = new IBond[2];
        System.arraycopy(ligandBonds, 0, arrayCopy, 0, 2);
        return arrayCopy;
	}

	/** {@inheritDoc} */
	@Override
	@TestMethod("testGetStereoBond")
	public IBond getStereoBond() {
		return this.stereoBond;
	}

	/** {@inheritDoc} */
	@Override
	@TestMethod("testGetStereo")
	public Conformation getStereo() {
		return this.stereo;
	}

    @TestMethod("testMap_Map_Map,testMap_Null_Map,testMap_Map_Map_NullElement,testMap_Map_Map_EmptyMapping")
    @Override
    public IDoubleBondStereochemistry map(Map<IAtom, IAtom> atoms, Map<IBond, IBond> bonds) {

        if(bonds == null)
            throw new IllegalArgumentException("null bond mapping provided");

        // map the double bond and the connected ligand bonds
        IBond   doubleBond  = stereoBond != null ? bonds.get(stereoBond) : null;
        IBond[] connected   = new IBond[ligandBonds.length];

        for(int i = 0; i < connected.length; i++){
            if(ligandBonds[i] != null)
                connected[i] = bonds.get(ligandBonds[i]);
        }

        return new DoubleBondStereochemistry(doubleBond, connected, stereo);
    }

	@Override
	public <T extends ICDKObject> T newInstance(Class<T> clazz,
			Object... params) throws IllegalArgumentException {
		// TODO Auto-generated method stub
		return null;
	}
}
