/*
 * ModelAveragingCenancestorTreeLikelihood.java
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

package dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.AFsequence;
import dr.evomodel.branchratemodel.CenancestorBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractGeneralFrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Samples a single crypt-division model that is used at every internal node.
 *
 * @author Diego Mallo
 */
public class ModelAveragingCenancestorTreeLikelihood extends CenancestorTreeLikelihood {

    private final Parameter divisionModel;
    private final ModelAveragingCenancestorLikelihoodCore modelAveragingLikelihoodCore;
    private final boolean storePartials;

    public ModelAveragingCenancestorTreeLikelihood(PatternList patternList,
                                                   TreeModel treeModel,
                                                   SiteModel siteModel,
                                                   CenancestorBranchRateModel branchRateModel,
                                                   TipStatesModel tipStatesModel,
                                                   Parameter cenancestorHeight,
                                                   Parameter cenancestorBranch,
                                                   AbstractGeneralFrequencyModel cenancestorFrequencyModel,
                                                   Parameter divisionModel,
                                                   boolean useAmbiguities,
                                                   boolean allowMissingTaxa,
                                                   boolean storePartials,
                                                   boolean forceJavaCore,
                                                   boolean forceRescaling,
                                                   boolean heightRules) {
        super(patternList, treeModel, siteModel, branchRateModel, tipStatesModel,
                cenancestorHeight, cenancestorBranch, cenancestorFrequencyModel, "default",
                useAmbiguities, allowMissingTaxa, storePartials, forceJavaCore, forceRescaling,
                heightRules, createLikelihoodCore(patternList, divisionModel),
                "Java cenancestor FlipFlop model averaging");

        ModelAveragingCenancestorLikelihoodCore.validateDivisionModel(divisionModel);
        divisionModel.addBounds(new Parameter.DefaultBounds(
                ModelAveragingCenancestorLikelihoodCore.MODEL_COUNT - 1,
                ModelAveragingCenancestorLikelihoodCore.IDENTITY, 1));
        this.divisionModel = divisionModel;
        this.modelAveragingLikelihoodCore =
                (ModelAveragingCenancestorLikelihoodCore) cenancestorlikelihoodCore;
        this.storePartials = storePartials;
        addVariable(divisionModel);
    }

    private static CenancestorLikelihoodCore createLikelihoodCore(PatternList patternList,
                                                                  Parameter divisionModel) {
        if (!(patternList.getDataType() instanceof AFsequence)) {
            throw new IllegalArgumentException(
                    "Crypt-division model averaging is only implemented for AFsequence data");
        }
        return new ModelAveragingCenancestorLikelihoodCore(
                patternList.getStateCount(), divisionModel);
    }

    public Parameter getDivisionModelParameter() {
        return divisionModel;
    }

    public int getDivisionModel() {
        return ModelAveragingCenancestorLikelihoodCore.getDivisionModelIndex(divisionModel);
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        if (variable == divisionModel) {
            modelAveragingLikelihoodCore.updateDivisionModel();
            updateAllNodes();
            fireModelChanged();
        } else {
            super.handleVariableChangedEvent(variable, index, type);
        }
    }

    @Override
    protected void storeState() {
        if (!storePartials) {
            modelAveragingLikelihoodCore.storeDivisionModelState();
        }
        super.storeState();
    }

    @Override
    protected void restoreState() {
        super.restoreState();
        if (!storePartials) {
            modelAveragingLikelihoodCore.restoreDivisionModelState();
        }
    }
}
