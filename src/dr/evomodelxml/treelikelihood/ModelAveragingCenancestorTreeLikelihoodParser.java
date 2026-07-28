/*
 * ModelAveragingCenancestorTreeLikelihoodParser.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.AFsequence;
import dr.evomodel.branchratemodel.CenancestorBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractGeneralFrequencyModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.ModelAveragingCenancestorTreeLikelihood;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses a cenancestor tree likelihood with a sampled crypt-division model.
 *
 * @author Diego Mallo
 */
public class ModelAveragingCenancestorTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String MODEL_AVERAGING_TREE_LIKELIHOOD =
            "modelAveragingCenancestorTreeLikelihood";
    public static final String DIVISION_MODEL = "divisionModel";

    public String getParserName() {
        return MODEL_AVERAGING_TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final boolean useAmbiguities = xo.getAttribute(
                CenancestorTreeLikelihoodParser.USE_AMBIGUITIES, false);
        final boolean allowMissingTaxa = xo.getAttribute(
                CenancestorTreeLikelihoodParser.ALLOW_MISSING_TAXA, false);
        final boolean storePartials = xo.getAttribute(
                CenancestorTreeLikelihoodParser.STORE_PARTIALS, true);
        boolean forceJavaCore = xo.getAttribute(
                CenancestorTreeLikelihoodParser.FORCE_JAVA_CORE, false);
        final boolean forceRescaling = xo.getAttribute(
                CenancestorTreeLikelihoodParser.FORCE_RESCALING, false);
        final boolean heightRules = xo.getAttribute(
                CenancestorTreeLikelihoodParser.HEIGHT_RULES, false);

        if (Boolean.valueOf(System.getProperty("java.only"))) {
            forceJavaCore = true;
        }

        final PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        if (!(patternList.getDataType() instanceof AFsequence)) {
            throw new XMLParseException(
                    MODEL_AVERAGING_TREE_LIKELIHOOD +
                            " is only implemented for AFsequence data");
        }

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
        final AbstractGeneralFrequencyModel cenancestorFrequencyModel =
                (AbstractGeneralFrequencyModel) xo.getChild(AbstractGeneralFrequencyModel.class);
        final CenancestorBranchRateModel branchRateModel =
                (CenancestorBranchRateModel) xo.getChild(CenancestorBranchRateModel.class);
        final Parameter cenancestor =
                (Parameter) xo.getElementFirstChild(
                        CenancestorTreeLikelihoodParser.CENANCESTOR_HEIGHT);
        final Parameter cenancestorBranch =
                (Parameter) xo.getElementFirstChild(
                        CenancestorTreeLikelihoodParser.CENANCESTOR_BRANCH);
        final Parameter divisionModel =
                (Parameter) xo.getElementFirstChild(DIVISION_MODEL);

        final TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);
        if (tipStatesModel != null && tipStatesModel.getPatternList() != null) {
            throw new XMLParseException(
                    "The same sequence error model cannot be used for multiple partitions");
        }
        if (tipStatesModel != null && tipStatesModel.getModelType() == TipStatesModel.Type.STATES) {
            throw new XMLParseException("The state emitting TipStateModel requires BEAGLE");
        }

        try {
            return new ModelAveragingCenancestorTreeLikelihood(
                    patternList, treeModel, siteModel, branchRateModel, tipStatesModel,
                    cenancestor, cenancestorBranch, cenancestorFrequencyModel, divisionModel,
                    useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore,
                    forceRescaling, heightRules);
        } catch (IllegalArgumentException exception) {
            throw new XMLParseException(exception.getMessage());
        }
    }

    public String getParserDescription() {
        return "A cenancestor tree likelihood that samples one crypt-division model " +
                "for all internal nodes.";
    }

    public Class getReturnType() {
        return ModelAveragingCenancestorTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.ALLOW_MISSING_TAXA, true),
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.STORE_PARTIALS, true),
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.FORCE_JAVA_CORE, true),
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.FORCE_RESCALING, true),
            AttributeRule.newBooleanRule(CenancestorTreeLikelihoodParser.HEIGHT_RULES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteModel.class),
            new ElementRule(CenancestorBranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            new ElementRule(FrequencyModel.class, true),
            new ElementRule(CenancestorTreeLikelihoodParser.CENANCESTOR_HEIGHT,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(CenancestorTreeLikelihoodParser.CENANCESTOR_BRANCH,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(DIVISION_MODEL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
