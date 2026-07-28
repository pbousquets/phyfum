/*
 * ModelAveragingCenancestorLikelihoodCoreTest.java
 *
 * This file is part of PHYFUM.
 */

package test.dr.evomodel.treelikelihood;

import dr.evomodel.treelikelihood.BuddingCenancestorLikelihoodCore;
import dr.evomodel.treelikelihood.CenancestorLikelihoodCore;
import dr.evomodel.treelikelihood.FissionCenancestorLikelihoodCore;
import dr.evomodel.treelikelihood.GeneralCenancestorLikelihoodCore;
import dr.evomodel.treelikelihood.ModelAveragingCenancestorLikelihoodCore;
import dr.evomodel.treelikelihood.SplitFissionCenancestorLikelihoodCore;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Tests the model selector against the four fixed pruning cores.
 */
public class ModelAveragingCenancestorLikelihoodCoreTest extends TestCase {

    private static final int STATE_COUNT = 6;
    private static final int PATTERN_COUNT = 2;
    private static final int PARENT_NODE = 2;
    private static final double TOLERANCE = 1.0E-14;

    public void testEachDivisionModelMatchesFixedCore() {
        final Parameter divisionModel = new Parameter.Default(0.0);

        final CenancestorLikelihoodCore[] fixedCores = new CenancestorLikelihoodCore[]{
                new GeneralCenancestorLikelihoodCore(STATE_COUNT),
                new BuddingCenancestorLikelihoodCore(STATE_COUNT),
                new FissionCenancestorLikelihoodCore(STATE_COUNT),
                new SplitFissionCenancestorLikelihoodCore(STATE_COUNT)
        };

        for (int model = 0; model < fixedCores.length; model++) {
            divisionModel.setParameterValue(0, model);
            final CenancestorLikelihoodCore averagingCore =
                    new ModelAveragingCenancestorLikelihoodCore(STATE_COUNT, divisionModel);

            assertEquals(calculatePartials(fixedCores[model]),
                    calculatePartials(averagingCore), TOLERANCE);
        }
    }

    public void testCoreStateRestoresAfterRejectedModelChange() {
        final Parameter divisionModel = new Parameter.Default(0.0);
        final CenancestorLikelihoodCore core =
                new ModelAveragingCenancestorLikelihoodCore(STATE_COUNT, divisionModel);

        initializeCore(core);
        core.calculatePartials(0, 1, PARENT_NODE);
        final double[] identityPartials = getPartials(core);

        divisionModel.storeParameterValues();
        core.storeState();

        divisionModel.setParameterValue(0, 1.0);
        core.setNodePartialsForUpdate(PARENT_NODE);
        core.calculatePartials(0, 1, PARENT_NODE);

        divisionModel.restoreParameterValues();
        core.restoreState();

        assertEquals(identityPartials, getPartials(core), TOLERANCE);
    }

    public void testInvalidDivisionModelIsRejected() {
        final Parameter divisionModel = new Parameter.Default(0.5);
        try {
            ModelAveragingCenancestorLikelihoodCore.getDivisionModelIndex(divisionModel);
            fail("Expected a non-integer divisionModel value to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("integer"));
        }

        divisionModel.setParameterValue(0, 4.0);
        try {
            ModelAveragingCenancestorLikelihoodCore.getDivisionModelIndex(divisionModel);
            fail("Expected an out-of-range divisionModel value to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("0 to 3"));
        }
    }

    private double[] calculatePartials(CenancestorLikelihoodCore core) {
        initializeCore(core);
        core.calculatePartials(0, 1, PARENT_NODE);
        return getPartials(core);
    }

    private void initializeCore(CenancestorLikelihoodCore core) {
        core.initialize(3, PATTERN_COUNT, 1, true);
        core.overridableInitialization();

        core.setNodePartials(0, new double[]{
                0.12, 0.20, 0.18, 0.10, 0.15, 0.25,
                0.30, 0.10, 0.20, 0.15, 0.05, 0.20
        });
        core.setNodePartials(1, new double[]{
                0.22, 0.08, 0.20, 0.16, 0.14, 0.20,
                0.10, 0.25, 0.15, 0.20, 0.18, 0.12
        });
        core.createNodePartials(PARENT_NODE);

        core.setNodeMatrix(0, 0, createTransitionMatrix(0.70));
        core.setNodeMatrix(1, 0, createTransitionMatrix(0.55));
    }

    private double[] createTransitionMatrix(double diagonal) {
        final double[] matrix = new double[STATE_COUNT * STATE_COUNT];
        final double offDiagonal = (1.0 - diagonal) / (STATE_COUNT - 1);
        for (int row = 0; row < STATE_COUNT; row++) {
            for (int column = 0; column < STATE_COUNT; column++) {
                matrix[row * STATE_COUNT + column] =
                        row == column ? diagonal : offDiagonal;
            }
        }
        return matrix;
    }

    private double[] getPartials(CenancestorLikelihoodCore core) {
        final double[] partials = new double[STATE_COUNT * PATTERN_COUNT];
        core.getPartials(PARENT_NODE, partials);
        return partials;
    }

    private void assertEquals(double[] expected, double[] actual, double tolerance) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("Partial " + i, expected[i], actual[i], tolerance);
        }
    }
}
