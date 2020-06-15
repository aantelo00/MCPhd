/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2008  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
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
package org.openscience.cdk.debug;

import java.util.Map;

import org.openscience.cdk.formula.MolecularFormula;
import org.openscience.cdk.interfaces.IIsotope;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

/**
 * @cdk.module datadebug
 */
public class DebugMolecularFormula extends MolecularFormula implements IMolecularFormula {

    private ILoggingTool logger =
        LoggingToolFactory.createLoggingTool(DebugMolecularFormula.class);

	public IMolecularFormula add(IMolecularFormula formula) {
	    logger.debug("Adding formula: ", formula);
	    return super.add(formula);
    }

	public IMolecularFormula addIsotope(IIsotope isotope) {
	    logger.debug("Adding isotope: ", isotope);
	    return super.addIsotope(isotope);
    }

	public IMolecularFormula addIsotope(IIsotope isotope, int count) {
		logger.debug("Adding 'count' isotopes: ", isotope);
	    return super.addIsotope(isotope, count);
    }

	public boolean contains(IIsotope isotope) {
		logger.debug("Contains isotope?: ", isotope);
	    return super.contains(isotope);
    }

	public Integer getCharge() {
		Integer charge = super.getCharge();
		logger.debug("Getting charge: ", charge);
	    return charge;
    }

	public int getIsotopeCount(IIsotope isotope) {
	    logger.debug("Getting isotope count for: ", isotope);
	    return super.getIsotopeCount(isotope);
    }

	public int getIsotopeCount() {
		logger.debug("Getting isotope count: ", super.getIsotopeCount());
	    return super.getIsotopeCount();
    }

	public Map<Object, Object> getProperties() {
	    logger.debug("Getting properties...");
      return super.getProperties();
    }

	public Object getProperty(Object description) {
	    logger.debug("Getting property: " + description);
	    return super.getProperty(description);
    }

	public Iterable<IIsotope> isotopes() {
		logger.debug("Getting isotope iterator..");
	    return super.isotopes();
    }

	public void removeAllIsotopes() {
	    logger.debug("Removing all isotopes...");
	    super.removeAllIsotopes();
    }

	public void removeIsotope(IIsotope isotope) {
	    logger.debug("Removing this isotope: ", isotope);
	    super.removeIsotope(isotope);
    }

	public void removeProperty(Object description) {
	    logger.debug("Removing property: " + description);
      super.removeProperty(description);
    }

	public void setCharge(Integer charge) {
		logger.debug("Setting the charge to: ", charge);
	    super.setCharge(charge);
    }

	public void setProperties(Map<Object, Object> properties) {
	    logger.debug("Setting new properties...");
	    super.setProperties(properties);
    }

	public void setProperty(Object description, Object property) {
	    logger.debug("Setting new property: " + description + " -> " + property);
	    super.setProperty(description, property);
    }

    public IChemObjectBuilder getBuilder() {
        return DebugChemObjectBuilder.getInstance();
    }

}
