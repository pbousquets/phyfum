/*
 * DivisionModelPriorTest.java
 *
 * This file is part of PHYFUM.
 */

package test.dr.evomodel.treelikelihood;

import dr.evomodel.treelikelihood.DivisionModelPrior;
import dr.inference.model.Parameter;
import dr.inference.operators.UniformIntegerOperator;
import junit.framework.TestCase;

/**
 * Tests the categorical prior on the crypt-division model.
 */
public class DivisionModelPriorTest extends TestCase {

    public void testEqualWeights() {
        final Parameter divisionModel = new Parameter.Default(0.0);
        final DivisionModelPrior prior = new DivisionModelPrior(
                divisionModel, new double[]{1.0, 1.0, 1.0, 1.0});

        for (int model = 0; model < 4; model++) {
            divisionModel.setParameterValue(0, model);
            assertEquals(Math.log(0.25), prior.getLogLikelihood(), 0.0);
        }
    }

    public void testUnequalWeights() {
        final Parameter divisionModel = new Parameter.Default(0.0);
        final DivisionModelPrior prior = new DivisionModelPrior(
                divisionModel, new double[]{1.0, 2.0, 3.0, 4.0});

        for (int model = 0; model < 4; model++) {
            divisionModel.setParameterValue(0, model);
            assertEquals(Math.log((model + 1.0) / 10.0),
                    prior.getLogLikelihood(), 0.0);
        }
    }

    public void testZeroWeightExcludesModel() {
        final Parameter divisionModel = new Parameter.Default(2.0);
        final DivisionModelPrior prior = new DivisionModelPrior(
                divisionModel, new double[]{1.0, 1.0, 0.0, 1.0});

        assertEquals(Double.NEGATIVE_INFINITY, prior.getLogLikelihood());
    }

    public void testUniformIntegerProposalHasZeroHastingsRatio() {
        final Parameter divisionModel = new Parameter.Default(0.0);
        divisionModel.addBounds(new Parameter.DefaultBounds(3.0, 0.0, 1));
        final UniformIntegerOperator operator =
                new UniformIntegerOperator(divisionModel, 0, 3, 1.0);

        for (int i = 0; i < 100; i++) {
            assertEquals(0.0, operator.doOperation(), 0.0);
            final double value = divisionModel.getParameterValue(0);
            assertEquals(Math.rint(value), value, 0.0);
            assertTrue(value >= 0.0 && value <= 3.0);
        }
    }

    public void testInvalidWeightsAreRejected() {
        assertInvalidWeights(new double[]{1.0, 1.0, 1.0});
        assertInvalidWeights(new double[]{1.0, -1.0, 1.0, 1.0});
        assertInvalidWeights(new double[]{0.0, 0.0, 0.0, 0.0});
        assertInvalidWeights(new double[]{1.0, Double.NaN, 1.0, 1.0});
        assertInvalidWeights(new double[]{
                Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE});
    }

    private void assertInvalidWeights(double[] weights) {
        try {
            new DivisionModelPrior(new Parameter.Default(0.0), weights);
            fail("Expected invalid weights to be rejected");
        } catch (IllegalArgumentException expected) {
            // Expected.
        }
    }
}
