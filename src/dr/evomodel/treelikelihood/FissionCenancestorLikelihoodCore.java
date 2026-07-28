/*
 * FissionCenancestorLikelihoodCore.java
 *
 * Diego Mallo
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
 *  PHYFUM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with PHYFUM; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treelikelihood;

import java.util.Arrays;
import dr.math.ThreeVHypergeometricDistribution;

/**
 * FissionCenancestorLikelihoodCore - An implementation of LikelihoodCore for small cell populations that reproduce by fission
 * This means that each cell in the niche first replicates to form a 2S niche, that than splits in two S populations at random
 *
 * @author Diego Mallo
 */

public class FissionCenancestorLikelihoodCore extends GeneralCenancestorLikelihoodCore {
    protected int S;
    protected int stateCount;
    protected double [][] pFiss; //[parentalState][fission] probability
    protected int [][][] combFiss;//[parentalState][daughter][fission] state

    //Vars for calculations so that we are not allocating memory all the time
    protected double [] pLA,pLB;

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public FissionCenancestorLikelihoodCore(int stateCount) {
        super(stateCount);
        this.stateCount = stateCount;
        this.S = (int) (Math.sqrt(2 * stateCount + 1.0/4.0)-3.0/2.0);
        pLA = new double[stateCount];
        pLB = new double[stateCount];
    }

    public void overridableInitialization() {
        generateFissionTables();
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3)
    {
        double sumAbud , sumBbud; //Likelihood sums of the valid budding options on the [A,B] branch
        int u = 0;
        int v = 0;

        for (int l = 0; l < matrixCount; l++) {

            for (int k = 0; k < patternCount; k++) {

                int w = l * matrixSize;

                Arrays.fill(pLA,0.0);
                Arrays.fill(pLB, 0.0);

                //Calculate the partials at the end of the branch of each daughter branch A,B
                for (int i = 0; i < stateCount; i++) { //i = from state
                    for (int j = 0; j < stateCount; j++) { // j = to state
                        pLA[i] += matrices1[w] * partials1[v + j];
                        pLB[i] += matrices2[w] * partials2[v + j];
                        w++;
                    }
                }

                //Calculate the partials for the paternal node, taking into account all possible combinations of starting states for the daughter lineages
                for (int iState = 0; iState < stateCount; iState++) {
                    partials3[u] = 0;
                    for (int iFission = 0; iFission < combFiss[iState][0].length; iFission++) {
                        partials3[u] += pLA[combFiss[iState][0][iFission]] * pLB[combFiss[iState][1][iFission]] * pFiss[iState][iFission];
                    }
                    u++;
                }
                v += stateCount;
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3, int[] matrixMap)
    {
        double sumAbud , sumBbud; //Likelihood sums of the valid budding options on the [A,B] branch

        int u = 0;
        int v = 0;

        for (int k = 0; k < patternCount; k++) {

            int w = matrixMap[k] * matrixSize;

            Arrays.fill(pLA,0.0);
            Arrays.fill(pLB, 0.0);

            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    pLA[i] += matrices1[w] * partials1[v + j];
                    pLB[i] += matrices2[w] * partials2[v + j];
                    w++;
                }
            }

            //Calculate the partials for the paternal node, taking into account all possible combinations of starting states for the daughter lineages
            for (int iState = 0; iState < stateCount; iState++) {
                partials3[u] = 0;
                for (int iFission = 0; iFission < combFiss[iState][0].length; iFission++) {
                    partials3[u] += pLA[combFiss[iState][0][iFission]] * pLB[combFiss[iState][1][iFission]] * pFiss[iState][iFission];
                }
                u++;
            }
            v += stateCount;
        }
    }

    //This would need to be ported to fission
    //We are assuming fission is not happening at the origin, so we don't change the implementation of these
    /*   *//**
     * Calculates partial likelihoods at a node when both children have partials.
     *//*
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials3)


    *//**
     * Calculates partial likelihoods at a node when both children have partials.
     *//*
    protected void calculatePartialsPruning(double[] partials1, double[] matrices1,
                                            double[] partials3, int[] matrixMap)*/

    /**
     * Pending to implement. These are not currently in use in PHYFUM but I should still implement them, just not right now
     */

    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */

    protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                      double[] partials2, double[] matrices2,
                                                      double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when both children have states.
     */
    protected void calculateStatesStatesPruning(int[] states1, double[] matrices1,
                                                int[] states2, double[] matrices2,
                                                double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(	int[] states1, double[] matrices1,
                                                      double[] partials2, double[] matrices2,
                                                      double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }
    /**
     * Calculates partial likelihoods at a node when the only child has states.
     */
    protected void calculateStatesPruning(int[] states1, double[] matrices1,
                                          double[] partials3)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }

    /**
     * Calculates partial likelihoods at a node when the only child has states.
     */
    protected void calculateStatesPruning(int[] states1, double[] matrices1,
                                          double[] partials3, int[] matrixMap)
    {
        throw new UnsupportedOperationException("This method is not implemented.");
    }


    /**
     * This method creates tables with the probability of each fission combination and what those are
     * For each state, we count the number of fully demethylated cells S-k-m, half-methylated k, and fully-methylated m
     * For S=2, possible states are for (k,m): 2,0,0 1,0,1 0,0,2 1,1,0, 0,2,0 0,1,1
     **/
    private void generateFissionTables(){
        ThreeVHypergeometricDistribution iStateCombs;
        int[][] iStateCounts;

        //Initializing array with -1 to generate errors is for some reason the code tries to use unavailable states
        int[][] varState = new int[S+1][S+1]; //[k][m]{state}
        for (int i = 0; i < S+1; i++){
            Arrays.fill(varState[i],-1);
        }

        this.pFiss = new double[this.stateCount][];//[parentalState][fission] probability
        this.combFiss = new int[this.stateCount][][];//[parentalState][daughter][fission] state


        //Using two passes because I need all states already defined in the second. This can probably be improved but
        //is run only once so it is not needed
        int iState = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k <= S; k++){
                if (k+m <=S) {
                    varState[k][m] = iState;
                    iState++;
                }
            }
        }
        //state loop that calculates the multivariate hypergeometric combinations and probabilities and stores them
        //in the proper arrays
        //helpers for readability
        int ka,ma,k0,m0,k1,m1;
        iState = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k+m <= S; k++){
                if (k+m <=S) {
                    ka = k * 2;
                    ma = m * 2;

                    //Generation of combinations
                    iStateCombs = new ThreeVHypergeometricDistribution((S - k - m) * 2, ka, ma, S);

                    //getters
                    this.pFiss[iState] = iStateCombs.getProbs();
                    iStateCounts = iStateCombs.getCounts();

                    //Memory allocation
                    this.combFiss[iState] = new int[2][];
                    this.combFiss[iState][0] = new int[iStateCounts.length];
                    this.combFiss[iState][1] = new int[iStateCounts.length];

                    // translating counts given by k and m into states
                    for (int iFiss = 0; iFiss < iStateCounts.length; iFiss++) {

                        k0 = iStateCounts[iFiss][1];
                        m0 = iStateCounts[iFiss][2];
                        k1 = ka - k0; //The cells that do not go to daughter 0 go to 1
                        m1 = ma - m0; //The cells that do not go to daughter 0 go to 1

                        combFiss[iState][0][iFiss] = varState[k0][m0];
                        combFiss[iState][1][iFiss] = varState[k1][m1];
                    }
                    iState++;
                }
            }
        }
    }
}

