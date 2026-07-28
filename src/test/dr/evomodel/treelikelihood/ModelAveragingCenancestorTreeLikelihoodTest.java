/*
 * ModelAveragingCenancestorTreeLikelihoodTest.java
 *
 * This file is part of PHYFUM.
 */

package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.StrictClockCenancestorBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.CenancestorTreeLikelihood;
import dr.evomodel.treelikelihood.ModelAveragingCenancestorTreeLikelihood;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Tests model switching in the complete cenancestor tree likelihood.
 */
public class ModelAveragingCenancestorTreeLikelihoodTest extends TestCase {

    private static final String[] FIXED_MODEL_NAMES =
            new String[]{"identity", "budding", "fission", "split"};

    public void testEachModelMatchesFixedTreeLikelihood() throws Exception {
        for (int model = 0; model < FIXED_MODEL_NAMES.length; model++) {
            final CenancestorTreeLikelihood fixedLikelihood =
                    createFixedLikelihood(FIXED_MODEL_NAMES[model]);
            final AveragingLikelihood averaging = createAveragingLikelihood(model);

            assertEquals("division model " + model,
                    fixedLikelihood.getLogLikelihood(),
                    averaging.likelihood.getLogLikelihood(), 1.0E-12);
        }
    }

    public void testLikelihoodRestoresAfterRejectedModelChange() throws Exception {
        final AveragingLikelihood averaging = createAveragingLikelihood(0);
        final double identityLogLikelihood = averaging.likelihood.getLogLikelihood();

        averaging.likelihood.storeModelState();
        averaging.divisionModel.setParameterValue(0, 1.0);
        final double buddingLogLikelihood = averaging.likelihood.getLogLikelihood();
        assertTrue(identityLogLikelihood != buddingLogLikelihood);

        averaging.likelihood.restoreModelState();

        assertEquals(0.0, averaging.divisionModel.getParameterValue(0), 0.0);
        assertEquals(identityLogLikelihood,
                averaging.likelihood.getLogLikelihood(), 0.0);
    }

    private CenancestorTreeLikelihood createFixedLikelihood(String divisionModel)
            throws Exception {
        final Components components = createComponents();
        return new CenancestorTreeLikelihood(
                components.patterns, components.treeModel, components.siteModel,
                components.branchRates, components.errorModel,
                new Parameter.Default("cenancestorHeight", 3.0, 0.0, 100.0),
                new Parameter.Default("cenancestorBranch", 1.0, 0.0, 100.0), null,
                divisionModel, false, false, true, true, false, false);
    }

    private AveragingLikelihood createAveragingLikelihood(int model) throws Exception {
        final Components components = createComponents();
        final Parameter divisionModel = new Parameter.Default((double) model);
        final ModelAveragingCenancestorTreeLikelihood likelihood =
                new ModelAveragingCenancestorTreeLikelihood(
                        components.patterns, components.treeModel, components.siteModel,
                        components.branchRates, components.errorModel,
                        new Parameter.Default("cenancestorHeight", 3.0, 0.0, 100.0),
                        new Parameter.Default("cenancestorBranch", 1.0, 0.0, 100.0), null,
                        divisionModel, false, false, true, true, false, false);
        return new AveragingLikelihood(likelihood, divisionModel);
    }

    private Components createComponents() throws Exception {
        final int stemCells = 2;
        final AFsequence dataType = new AFsequence(6);
        final Taxa taxa = new Taxa();
        for (int i = 0; i < 3; i++) {
            taxa.addTaxon(new Taxon("C" + i));
        }

        final String[] sequences = new String[]{
                "0.1, 0.4, 0.14, 0.12",
                "0.2, 0.3, 0.13, 0.42",
                "0.6, 0.36, 0.1, 0.62"
        };
        final List<int[]> states = new ArrayList<int[]>();
        for (String sequence : sequences) {
            states.add(new AFsequence(sequence).getSequence());
        }

        final Patterns patterns = new Patterns(dataType, taxa);
        for (int site = 0; site < states.get(0).length; site++) {
            final int[] pattern = new int[states.size()];
            for (int taxon = 0; taxon < states.size(); taxon++) {
                pattern[taxon] = states.get(taxon)[site];
            }
            patterns.addPattern(pattern);
        }

        final Tree tree = new NewickImporter(
                "(C0:2.0,(C1:1.0,C2:1.0):1.0);").importTree(null);
        final TreeModel treeModel = new TreeModel(tree);

        final Parameter stemCellParameter = new Parameter.Default((double) stemCells);
        final TipStatesModel errorModel = new DeterministicTipPartials();

        final double[] frequencies = new double[]{
                1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0,
                1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0
        };
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, frequencies);
        final FlipFlopModel substitutionModel = new FlipFlopModel(
                "test", dataType, stemCellParameter,
                new Parameter.Default(0.05), new Parameter.Default(0.95),
                new Parameter.Default(0.05), false, true, frequencyModel);
        final GammaSiteModel siteModel = new GammaSiteModel(substitutionModel);
        final StrictClockCenancestorBranchRates branchRates =
                new StrictClockCenancestorBranchRates(new Parameter.Default(1.0));

        return new Components(patterns, treeModel, siteModel, errorModel, branchRates);
    }

    private static final class Components {
        private final Patterns patterns;
        private final TreeModel treeModel;
        private final GammaSiteModel siteModel;
        private final TipStatesModel errorModel;
        private final StrictClockCenancestorBranchRates branchRates;

        private Components(Patterns patterns, TreeModel treeModel, GammaSiteModel siteModel,
                           TipStatesModel errorModel,
                           StrictClockCenancestorBranchRates branchRates) {
            this.patterns = patterns;
            this.treeModel = treeModel;
            this.siteModel = siteModel;
            this.errorModel = errorModel;
            this.branchRates = branchRates;
        }
    }

    private static final class DeterministicTipPartials extends TipStatesModel {

        private DeterministicTipPartials() {
            super("deterministicTipPartials", null, null);
        }

        public Type getModelType() {
            return Type.PARTIALS;
        }

        public void getTipPartials(int nodeIndex, double[] partials) {
            for (int i = 0; i < partials.length; i++) {
                partials[i] = 0.1 + ((nodeIndex + i) % 7) * 0.1;
            }
        }

        public void getTipStates(int nodeIndex, int[] tipStates) {
            throw new IllegalArgumentException("This model emits only tip partials");
        }

        protected void taxaChanged() {
        }

        protected void handleModelChangedEvent(Model model, Object object, int index) {
            fireModelChanged();
        }

        protected void handleVariableChangedEvent(Variable variable, int index,
                                                  Parameter.ChangeType type) {
            fireModelChanged();
        }

        protected void storeState() {
        }

        protected void restoreState() {
        }

        protected void acceptState() {
        }
    }

    private static final class AveragingLikelihood {
        private final ModelAveragingCenancestorTreeLikelihood likelihood;
        private final Parameter divisionModel;

        private AveragingLikelihood(ModelAveragingCenancestorTreeLikelihood likelihood,
                                    Parameter divisionModel) {
            this.likelihood = likelihood;
            this.divisionModel = divisionModel;
        }
    }
}
