/*
 * DivisionModelPrior.java
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

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * A categorical prior on the crypt-division model indicator.
 *
 * The supplied values are relative weights and are normalized internally.
 *
 * @author Diego Mallo
 */
public class DivisionModelPrior extends AbstractModelLikelihood {

    public static final String DIVISION_MODEL_PRIOR = "divisionModelPrior";

    private final Parameter divisionModel;
    private final double[] logProbabilities;

    public DivisionModelPrior(Parameter divisionModel, double[] weights) {
        super(DIVISION_MODEL_PRIOR);

        if (divisionModel == null) {
            throw new IllegalArgumentException("divisionModel parameter cannot be null");
        }
        ModelAveragingCenancestorLikelihoodCore.validateDivisionModel(divisionModel);

        if (weights == null ||
                weights.length != ModelAveragingCenancestorLikelihoodCore.MODEL_COUNT) {
            throw new IllegalArgumentException("divisionModel prior requires four weights");
        }

        double sum = 0.0;
        for (double weight : weights) {
            if (Double.isNaN(weight) || Double.isInfinite(weight) || weight < 0.0) {
                throw new IllegalArgumentException(
                        "divisionModel prior weights must be finite and non-negative");
            }
            sum += weight;
        }
        if (sum <= 0.0 || Double.isInfinite(sum)) {
            throw new IllegalArgumentException(
                    "divisionModel prior weights must have a finite positive sum");
        }

        logProbabilities = new double[weights.length];
        for (int i = 0; i < weights.length; i++) {
            logProbabilities[i] = weights[i] == 0.0
                    ? Double.NEGATIVE_INFINITY
                    : Math.log(weights[i] / sum);
        }

        this.divisionModel = divisionModel;
        addVariable(divisionModel);
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return logProbabilities[
                ModelAveragingCenancestorLikelihoodCore.getDivisionModelIndex(divisionModel)];
    }

    public void makeDirty() {
        // The prior is inexpensive and is calculated directly on every request.
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // No submodels.
    }

    protected void handleVariableChangedEvent(Variable variable, int index,
                                              Parameter.ChangeType type) {
        // getLogLikelihood reads the current value directly.
    }

    protected void storeState() {
        // No derived state.
    }

    protected void restoreState() {
        // No derived state.
    }

    protected void acceptState() {
        // No derived state.
    }
}
