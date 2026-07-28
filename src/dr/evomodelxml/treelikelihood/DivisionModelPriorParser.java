/*
 * DivisionModelPriorParser.java
 *
 * This file is part of PHYFUM.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * PHYFUM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * PHYFUM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

package dr.evomodelxml.treelikelihood;

import dr.evomodel.treelikelihood.DivisionModelPrior;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses a categorical prior on the crypt-division model.
 *
 * @author Diego Mallo
 */
public class DivisionModelPriorParser extends AbstractXMLObjectParser {

    public static final String WEIGHTS = "weights";

    public String getParserName() {
        return DivisionModelPrior.DIVISION_MODEL_PRIOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final Parameter divisionModel = (Parameter) xo.getChild(Parameter.class);
        final double[] weights = xo.hasAttribute(WEIGHTS)
                ? xo.getDoubleArrayAttribute(WEIGHTS)
                : new double[]{1.0, 1.0, 1.0, 1.0};

        try {
            return new DivisionModelPrior(divisionModel, weights);
        } catch (IllegalArgumentException exception) {
            throw new XMLParseException(exception.getMessage());
        }
    }

    public String getParserDescription() {
        return "A categorical prior on identity, budding, fission, and split-fission " +
                "crypt-division models.";
    }

    public Class getReturnType() {
        return DivisionModelPrior.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleArrayRule(WEIGHTS, true),
            new ElementRule(Parameter.class)
    };
}
