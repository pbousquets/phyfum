/*
 * ModelAveragingCenancestorLikelihoodCore.java
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

import dr.inference.model.Parameter;

/**
 * A cenancestor likelihood core that selects one crypt-division pruning
 * strategy for all internal nodes. The partials and matrices are stored only
 * once in this core; the delegate cores provide the pruning calculations.
 *
 * @author Diego Mallo
 */
public class ModelAveragingCenancestorLikelihoodCore extends GeneralCenancestorLikelihoodCore {

    public static final int IDENTITY = 0;
    public static final int BUDDING = 1;
    public static final int FISSION = 2;
    public static final int SPLIT_FISSION = 3;
    public static final int MODEL_COUNT = 4;

    private final Parameter divisionModel;
    private final GeneralCenancestorLikelihoodCore[] divisionCores;

    public ModelAveragingCenancestorLikelihoodCore(int stateCount, Parameter divisionModel) {
        super(stateCount);

        if (divisionModel == null) {
            throw new IllegalArgumentException("divisionModel parameter cannot be null");
        }

        this.divisionModel = divisionModel;
        divisionCores = new GeneralCenancestorLikelihoodCore[]{
                new GeneralCenancestorLikelihoodCore(stateCount),
                new BuddingCenancestorLikelihoodCore(stateCount),
                new FissionCenancestorLikelihoodCore(stateCount),
                new SplitFissionCenancestorLikelihoodCore(stateCount)
        };
    }

    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories) {
        super.initialize(nodeCount, patternCount, matrixCount, integrateCategories);

        // The delegates use these dimensions in their pruning loops, but do not
        // need their own partial or matrix buffers.
        for (GeneralCenancestorLikelihoodCore core : divisionCores) {
            core.nodeCount = nodeCount;
            core.patternCount = patternCount;
            core.matrixCount = matrixCount;
            core.integrateCategories = integrateCategories;
            core.partialsSize = partialsSize;
            core.matrixSize = matrixSize;
        }
    }

    @Override
    public void overridableInitialization() {
        for (GeneralCenancestorLikelihoodCore core : divisionCores) {
            core.overridableInitialization();
        }
    }

    public static int getDivisionModelIndex(Parameter divisionModel) {
        if (divisionModel.getDimension() != 1) {
            throw new IllegalArgumentException("divisionModel must have dimension 1");
        }

        final double value = divisionModel.getParameterValue(0);
        final int model = (int) value;

        if (value != model || model < IDENTITY || model >= MODEL_COUNT) {
            throw new IllegalArgumentException(
                    "divisionModel must be an integer from 0 to " + (MODEL_COUNT - 1) + ": " + value);
        }
        return model;
    }

    private GeneralCenancestorLikelihoodCore getDivisionCore() {
        return divisionCores[getDivisionModelIndex(divisionModel)];
    }

    @Override
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3) {
        getDivisionCore().calculateStatesStatesPruning(
                states1, matrices1, states2, matrices2, partials3);
    }

    @Override
    protected void calculateStatesPartialsPruning(int[] states1, double[] matrices1,
                                                  double[] partials2, double[] matrices2,
                                                  double[] partials3) {
        getDivisionCore().calculateStatesPartialsPruning(
                states1, matrices1, partials2, matrices2, partials3);
    }

    @Override
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3) {
        getDivisionCore().calculatePartialsPartialsPruning(
                partials1, matrices1, partials2, matrices2, partials3);
    }

    @Override
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3, int[] matrixMap) {
        getDivisionCore().calculateStatesStatesPruning(
                states1, matrices1, states2, matrices2, partials3, matrixMap);
    }

    @Override
    protected void calculateStatesPartialsPruning(int[] states1, double[] matrices1,
                                                  double[] partials2, double[] matrices2,
                                                  double[] partials3, int[] matrixMap) {
        getDivisionCore().calculateStatesPartialsPruning(
                states1, matrices1, partials2, matrices2, partials3, matrixMap);
    }

    @Override
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3, int[] matrixMap) {
        getDivisionCore().calculatePartialsPartialsPruning(
                partials1, matrices1, partials2, matrices2, partials3, matrixMap);
    }

    @Override
    protected void calculateStatesPruning(int[] states1, double[] matrices1, double[] partials3) {
        getDivisionCore().calculateStatesPruning(states1, matrices1, partials3);
    }

    @Override
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                            double[] partials3) {
        getDivisionCore().calculatePartialsPruning(partials1, matrices1, partials3);
    }

    @Override
    protected void calculateStatesPruning(int[] states1, double[] matrices1,
                                          double[] partials3, int[] matrixMap) {
        getDivisionCore().calculateStatesPruning(states1, matrices1, partials3, matrixMap);
    }

    @Override
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                            double[] partials3, int[] matrixMap) {
        getDivisionCore().calculatePartialsPruning(partials1, matrices1, partials3, matrixMap);
    }

    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions,
                                              double[] outPartials) {
        getDivisionCore().calculateIntegratePartials(inPartials, proportions, outPartials);
    }
}
