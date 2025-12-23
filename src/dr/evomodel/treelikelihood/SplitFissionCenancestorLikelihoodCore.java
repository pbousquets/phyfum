/*
 * SplitFissionCenancestorLikelihoodCore.java
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
import jdk.jshell.spi.ExecutionControl;

import static java.lang.Math.ceil;

/**
 * SplitFissionCenancestorLikelihoodCore - An implementation of LikelihoodCore for small cell populations that reproduce
 * by splitting
 * If S is even, the niche is split at random in two populations of size S/2 and then grow back to S by instantaneous duplication of
 * all cells.
 * Not implemented for odd S.
 *
 * @author Diego Mallo
 */

public class SplitFissionCenancestorLikelihoodCore extends FissionCenancestorLikelihoodCore {

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public SplitFissionCenancestorLikelihoodCore(int stateCount) {
        super(stateCount);
    }

    public void overridableInitialization() {
        if(this.S%2==0) {
            generateEvenSplitFissionTables();
        } else {
            //generateOddSplitFissionTables();
        }
    }

    /**
     * This method creates tables with the probability of each split combination and what those are
     * For each state, we count the number of fully demethylated cells S-k-m, half-methylated k, and fully-methylated m
     * For S=2, possible states are for (k,m): 2,0,0 1,0,1 0,0,2 1,1,0, 0,2,0 0,1,1
     **/
    private void generateEvenSplitFissionTables(){
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
        int k0,m0,k1,m1;
        iState = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k+m <= S; k++){
                if (k+m <=S) {

                    //Generation of combinations
                    iStateCombs = new ThreeVHypergeometricDistribution(S - k - m, k, m, S/2);

                    //getters
                    this.pFiss[iState] = iStateCombs.getProbs();
                    iStateCounts = iStateCombs.getCounts();

                    //Memory allocation
                    this.combFiss[iState] = new int[2][];
                    this.combFiss[iState][0] = new int[iStateCounts.length];
                    this.combFiss[iState][1] = new int[iStateCounts.length];

                    // translating counts given by k and m into states
                    for (int iFiss = 0; iFiss < iStateCounts.length; iFiss++) {

                        k0 = iStateCounts[iFiss][1] * 2;
                        m0 = iStateCounts[iFiss][2] * 2;
                        k1 = k * 2 - k0; //The cells that do not go to daughter 0 go to 1
                        m1 = m * 2 - m0; //The cells that do not go to daughter 0 go to 1

                        combFiss[iState][0][iFiss] = varState[k0][m0];
                        combFiss[iState][1][iFiss] = varState[k1][m1];
                    }
                    iState++;
                }
            }
        }


    }

    /**
     * This method creates tables with the probability of each split fission combination and what those are
     * For each state, we count the number of fully demethylated cells S-k-m, half-methylated k, and fully-methylated m
     * For odd S this is much more complicated, since the daughter that takes ceil(S/2) needs to lose 1 cell at random
     * and the other gain one
     * For S=2, possible states are for (k,m): 2,0,0 1,0,1 0,0,2 1,1,0, 0,2,0 0,1,1
     **/
    //TODO splitting odd and even for development but I may merge them (or not) in the future
    private void generateOddSplitFissionTables(){
        throw new UnsupportedOperationException("SplitFission is not yet implemented for odd S");
    }
}

