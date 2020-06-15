/*  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 2006-2007  Egon Willighagen <egonw@users.sf.net>
 *
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
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
package org.openscience.cdk.nonotify;

import org.openscience.cdk.Bond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectListener;
import org.openscience.cdk.interfaces.IChemObjectBuilder;

/**
 * @cdk.module nonotify
 * @cdk.githash
 * @deprecated    Use the {@link org.openscience.cdk.silent.Bond} instead.
 */
public class NNBond extends Bond {

	private static final long serialVersionUID = -6286559199533034379L;

	public NNBond() {
		this(null, null, null, IBond.Stereo.NONE);
		this.atomCount = 0;
	}

	public NNBond(IAtom atom1, IAtom atom2) {
		this(atom1, atom2, IBond.Order.SINGLE, IBond.Stereo.NONE);
	}

	public NNBond(IAtom atom1, IAtom atom2, IBond.Order order) {
		this(atom1, atom2, order, IBond.Stereo.NONE);
	}

	public NNBond(IAtom atom1, IAtom atom2, IBond.Order order,
			      IBond.Stereo stereo) {
		super(atom1, atom2, order, stereo);
		setNotification(false);
	}

    public NNBond(IAtom[] atoms) {
        super(atoms);
        setNotification(false);
    }

    public NNBond(IAtom[] atoms, IBond.Order order) {
        super(atoms, order);
        setNotification(false);
    }

    public IChemObjectBuilder getBuilder() {
		return NoNotificationChemObjectBuilder.getInstance();
	}
	
	public void addListener(IChemObjectListener col) {
		// Ignore this: we do not listen anyway
	}
}

