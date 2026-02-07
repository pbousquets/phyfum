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

import java.util.*;

import dr.math.ThreeVHypergeometricDistribution;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;

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
            generateOddSplitFissionTables();
        }
    }

    /**
     * This method creates tables with the probability of each split combination and what those are
     * For each state, we count the number of fully demethylated cells S-k-m, half-methylated k, and fully-methylated m
     * For S=2, possible states are for (k,m): 2,0,0 1,0,1 0,0,2 1,1,0, 0,2,0 0,1,1
     **/
    private void generateEvenSplitFissionTables(){
        ThreeVHypergeometricDistribution iStateCombs;
        int[][] iStateCounts; //[combination]{d,k,m}

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
                    this.combFiss[iState] = new int[2][iStateCounts.length];

                    // translating counts given by k and m into states
                    for (int iFiss = 0; iFiss < iStateCounts.length; iFiss++) {

                        k0 = iStateCounts[iFiss][1] * 2;
                        m0 = iStateCounts[iFiss][2] * 2;
                        k1 = k * 2 - k0; //The cells that do not go to daughter 0 go to 1
                        m1 = m * 2 - m0; //The cells that do not go to daughter 0 go to 1

                        this.combFiss[iState][0][iFiss] = varState[k0][m0];
                        this.combFiss[iState][1][iFiss] = varState[k1][m1];
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
        ThreeVHypergeometricDistribution iStateCombs;
        SplitFissionCombinations iStateSplitFissionCombs;
        int[][][] iStateCounts; //[combination][daughter]{d,k,m}

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
        iState = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k+m <= S; k++){
                if (k+m <=S) {

                    //Generation of combinations
                    //Splitting combinations
                    iStateCombs = new ThreeVHypergeometricDistribution(S - k - m, k, m, (int) ceil(S/2.0)); //Large daughter

                    //getters
                    iStateSplitFissionCombs = new SplitFissionCombinations(iStateCombs);
                    iStateSplitFissionCombs.simplify();
                    iStateCounts = iStateSplitFissionCombs.getCounts();
                    this.pFiss[iState] = iStateSplitFissionCombs.getProbs();

                    //Memory allocation
                    this.combFiss[iState] = new int[2][iStateCounts.length];

                    // translating counts given by k and m into states
                    for (int iFiss = 0; iFiss < iStateCounts.length; iFiss++) {
                        this.combFiss[iState][0][iFiss] = varState[iStateCounts[iFiss][0][1]][iStateCounts[iFiss][0][2]];
                        this.combFiss[iState][1][iFiss] = varState[iStateCounts[iFiss][1][1]][iStateCounts[iFiss][1][2]];
                    }
                    iState++;
                }
            }
        }
    }

    private static class SplitFissionCombinations {

        int [][][] counts;//[combination][daughter]{d,k,m}
        double [] logProbs;
        double [] splitLogProbs;
        static final double precission = 1e-16;

        public SplitFissionCombinations(ThreeVHypergeometricDistribution iDistLarge) {
            this(iDistLarge, true);
        }

        public SplitFissionCombinations(ThreeVHypergeometricDistribution iDistLarge, boolean simplify) {
            int [][] iCombs = iDistLarge.getCounts();
            double [] iLps = iDistLarge.getLogProbs();
            List<int[]> combsBuilderA = new ArrayList<>();
            List<int[]> combsBuilderB = new ArrayList<>();
            List<Double> lpsBuilder = new ArrayList<>();
            List<Double> splitLpsBuilder = new ArrayList<>();

            int [] Ka = iDistLarge.getCollection();

            //Recalculate S
            int S = 0;
            for (int count : Ka) {
                S += count;
            }
            int St0a = (int) ceil(S/2.0)*2; //Number of cells in the large stem cell niche after split and duplication, before killing of 1
            int St0b = (int) floor(S/2.0)*2; //Number of cells in the small stem cell niche after split and duplication, secondary duplication of 1

            //Helpers
            int [] countsB = new int[3];
            int [] modCountsA = new int[3];
            int [] modCountsB = new int[3];
            int [] countsA = new int[3];
            double lp;

            //For each 3DHypergeometric distribution (split options) we calculate all possible alternative adjustments
            //(i.e., gain of 1 cell in the small subpopulation, loss of 1 in the big
            for (int iComb = 0; iComb < iLps.length; iComb++){
                System.arraycopy(iCombs[iComb],0,countsA,0,iCombs[iComb].length);

                for (int iCount = 0; iCount < countsA.length; iCount++) {
                    countsB[iCount] = (Ka[iCount] - countsA[iCount]) * 2; //Number of d,k,m cells in the small subpopulation after doubling
                    countsA[iCount] = countsA[iCount] * 2; //Number of d,k,m cells in the large subpopulation after doubling
                }

                for (int iCountA = 0; iCountA < 3; iCountA++) {
                    if(countsA[iCountA]==0) {
                        continue;
                    } else {
                        modCountsA = new int[3];
                        System.arraycopy(countsA,0,modCountsA,0,countsA.length);
                        modCountsA[iCountA] -= 1; //Remove one cell
                    }
                    for (int iCountB = 0; iCountB < 3; iCountB++) {
                        if(countsB[iCountB]==0){
                            continue;
                        } else {
                            modCountsB = new int[3];
                            System.arraycopy(countsB,0,modCountsB,0,countsB.length);
                            modCountsB[iCountB] += 1; //Add one cell
                        }

                        //Add the calculated combinations
                        combsBuilderA.add(modCountsA);
                        combsBuilderB.add(modCountsB);
                        lp = Math.log(0.5)+iLps[iComb]+Math.log(countsA[iCountA]/(double)St0a)+Math.log(countsB[iCountB]/(double)St0b); //1/2 of the two possibilities of assigning big/small S to each daughter lineage * split prob * rearrangement prob A * rearrangement prob B
                        lpsBuilder.add(lp);
                        splitLpsBuilder.add(iLps[iComb]);

                        //Add the same combination reversed, so that A is the small and B is the big
                        combsBuilderA.add(modCountsB);
                        combsBuilderB.add(modCountsA);
                        lpsBuilder.add(lp);
                        splitLpsBuilder.add(iLps[iComb]);
                    }
                }
            }

            int n = combsBuilderA.size();
            this.counts = new int[n][2][3];
            this.logProbs = new double[n];
            this.splitLogProbs = new double[n];

            for (int iComb = 0; iComb < n; iComb++) {
                System.arraycopy(combsBuilderA.get(iComb), 0, this.counts[iComb][0], 0, 3);
                System.arraycopy(combsBuilderB.get(iComb), 0, this.counts[iComb][1], 0, 3);
                this.logProbs[iComb] = lpsBuilder.get(iComb);
                this.splitLogProbs[iComb] = splitLpsBuilder.get(iComb);
            }

            if (simplify) {
                simplify();
            }
        }

        /** deep copy to avoid external mutation */
        public int[][][] getCounts() { return deepCopyCounts(counts); }

        /** return probabilities computed from logProbs via exp(logP) */
        public double[] getProbs() {
            double[] p = new double[logProbs.length];
            for (int i = 0; i < logProbs.length; i++) {
                p[i] = Math.exp(logProbs[i]);
            }
            // renormalize to remove tiny rounding drift so sum==1 within machine precision
            renormalizeIfNeeded(p);
            return p;
        }

        /** natural-log probabilities */
        public double[] getLogProbs() { return logProbs.clone(); }

        /** natural-log probabilities of the split (ThreeVHypergeometricDistribution), may be null after simplify */
        public double[] getSplitLogProbs() { return splitLogProbs == null ? null : splitLogProbs.clone(); }

        /** probabilities of the split (ThreeVHypergeometricDistribution), may be null after simplify */
        public double[] getSplitProbs() {
            if (splitLogProbs == null) return null;
            double[] p = new double[splitLogProbs.length];
            for (int i = 0; i < splitLogProbs.length; i++) {
                p[i] = Math.exp(splitLogProbs[i]);
            }
            return p;
        }

        private static int[][][] deepCopyCounts(int[][][] src) {
            int[][][] dst = new int[src.length][2][3];
            for (int i = 0; i < src.length; i++) {
                System.arraycopy(src[i][0], 0, dst[i][0], 0, 3);
                System.arraycopy(src[i][1], 0, dst[i][1], 0, 3);
            }
            return dst;
        }

        private static void renormalizeIfNeeded(double[] probs) {
            double sum = 0.0;
            for (double v : probs) sum += v;
            double err = Math.abs(sum - 1.0);
            if (err > precission && sum > 0.0) {
                for (int i = 0; i < probs.length; i++) probs[i] /= sum;
            }
        }

        public void simplify(){
            Map<CombKey, List<Double>> uniqueCounts = new HashMap<>();
            CombKey countKey;

            for (int iCount =0; iCount < this.counts.length; iCount++) {
                countKey = new CombKey(this.counts[iCount][0],this.counts[iCount][1],this.logProbs[iCount]);
                uniqueCounts.computeIfAbsent(countKey, k -> new ArrayList<>()).add(countKey.getLp());
            }

            int[][][] finalCombs = new int[uniqueCounts.size()][2][3];
            double[] finalLps = new double[uniqueCounts.size()];

            int iComb=0;
            for (Map.Entry<CombKey, List<Double>> e : uniqueCounts.entrySet()) {
                finalCombs[iComb][0] = e.getKey().getA();
                finalCombs[iComb][1] = e.getKey().getB();
                finalLps[iComb] = logSumList(e.getValue());
                iComb++;
            }

            this.counts = finalCombs;
            this.logProbs = finalLps;
            this.splitLogProbs = null;

        }

        //Helper class to hash int arrays and compare them by value for simplify
        private static class CombKey {
            private final int [] cArray;
            private final int Alength;
            private final int Blength;
            private final double lp;

            public CombKey(int[] arrayA, int[] arrayB, double lp){
                this.cArray = new int[arrayA.length+arrayB.length];
                System.arraycopy(arrayA,0,this.cArray,0,arrayA.length);
                System.arraycopy(arrayB,0,this.cArray,arrayA.length,arrayB.length);
                this.Alength = arrayA.length;
                this.Blength = arrayB.length;
                this.lp = lp;
            }

            public int [] getArray() {
                int [] retArray = new int [this.cArray.length];
                System.arraycopy(this.cArray,0,retArray,0,this.cArray.length);
                return(retArray);
            }

            public int [] getA() {
                int [] A = new int[this.Alength];
                System.arraycopy(this.cArray,0,A,0,this.Alength);
                return(A);
            }

            public int [] getB() {
                int [] B = new int [this.Blength];
                System.arraycopy(this.cArray,this.Alength,B,0,this.Blength);
                return(B);
            }

            public double getLp(){
                return(this.lp);
            }

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                if (!(o instanceof CombKey)) return false;
                CombKey other = (CombKey) o;
                return Arrays.equals(this.cArray, other.getArray());
            }

            @Override
            public int hashCode() {
                return Arrays.hashCode(this.cArray);
            }
        }

        static double logAdd(double a, double b) {
            // returns log(e^a + e^b) robustly
            if (Double.isInfinite(a) && a < 0) return b; // -Inf + b => b
            if (Double.isInfinite(b) && b < 0) return a;
            if (a < b) { double t = a; a = b; b = t; } // ensure a >= b
            return a + Math.log1p(Math.exp(b - a));
        }

        static double logSumList(List<Double> logps) {
            double acc = Double.NEGATIVE_INFINITY; // log(0)
            for (double lp : logps) {
                acc = logAdd(acc, lp);
            }
            return acc; // log(sum of probabilities)
        }
    }

    public static void printOddSplitFissionTables(int S, Integer filterK, Integer filterM) {
        int stateCount = (S + 1) * (S + 2) / 2;
        SplitFissionCenancestorLikelihoodCore core = new SplitFissionCenancestorLikelihoodCore(stateCount);
        if (core.S != S) {
            throw new IllegalArgumentException("Computed S mismatch. Expected " + S + " but got " + core.S);
        }
        core.overridableInitialization();

        int[] kByState = new int[stateCount];
        int[] mByState = new int[stateCount];
        int[][] varState = new int[S + 1][S + 1];
        for (int i = 0; i < S + 1; i++) {
            Arrays.fill(varState[i], -1);
        }
        int iState = 0;
        for (int m = 0; m <= S; m++) {
            for (int k = 0; k + m <= S; k++) {
                kByState[iState] = k;
                mByState[iState] = m;
                varState[k][m] = iState;
                iState++;
            }
        }

        Integer filterState = null;
        if (filterK != null && filterM != null) {
            if (filterK < 0 || filterM < 0 || filterK + filterM > S) {
                throw new IllegalArgumentException("Invalid (k,m) filter: k=" + filterK + " m=" + filterM);
            }
            int idx = 0;
            for (int m = 0; m <= S; m++) {
                for (int k = 0; k + m <= S; k++) {
                    if (k == filterK && m == filterM) {
                        filterState = idx;
                        break;
                    }
                    idx++;
                }
                if (filterState != null) break;
            }
        }

        for (int s = 0; s < stateCount; s++) {
            if (filterState != null && s != filterState) continue;
            int k = kByState[s];
            int m = mByState[s];
            int d = S - k - m;
            if (S % 2 == 0) {
                ThreeVHypergeometricDistribution iStateCombs = new ThreeVHypergeometricDistribution(d, k, m, S / 2);
                int[][] counts = iStateCombs.getCounts();
                double[] splitPs = iStateCombs.getProbs();
                double sum = 0.0;
                System.out.println("State " + s + " (d,k,m)=(" + d + "," + k + "," + m + ")");
                for (int iFiss = 0; iFiss < counts.length; iFiss++) {
                    int k0 = counts[iFiss][1] * 2;
                    int m0 = counts[iFiss][2] * 2;
                    int k1 = k * 2 - k0;
                    int m1 = m * 2 - m0;
                    int s0 = (k0 >= 0 && m0 >= 0 && k0 + m0 <= S) ? varState[k0][m0] : -1;
                    int s1 = (k1 >= 0 && m1 >= 0 && k1 + m1 <= S) ? varState[k1][m1] : -1;
                    double p = splitPs[iFiss];
                    sum += p;
                    String a;
                    if (s0 >= 0 && s0 < stateCount) {
                        int k0s = kByState[s0];
                        int m0s = mByState[s0];
                        int d0 = S - k0s - m0s;
                        a = "A=" + s0 + " (" + d0 + "," + k0 + "," + m0 + ")";
                    } else {
                        a = "A=INVALID(" + s0 + ")";
                    }
                    String b;
                    if (s1 >= 0 && s1 < stateCount) {
                        int k1s = kByState[s1];
                        int m1s = mByState[s1];
                        int d1 = S - k1s - m1s;
                        b = "B=" + s1 + " (" + d1 + "," + k1 + "," + m1 + ")";
                    } else {
                        b = "B=INVALID(" + s1 + ")";
                    }
                    System.out.printf("  p=%.12g  pSplit=%.12g  %s  %s%n", p, p, a, b);
                }
                System.out.printf("  sum(p)=%.15g%n", sum);
                System.out.println();
            } else {
                ThreeVHypergeometricDistribution iStateCombs = new ThreeVHypergeometricDistribution(d, k, m, (int) ceil(S / 2.0));
                SplitFissionCombinations rawCombs = new SplitFissionCombinations(iStateCombs, false);
                int[][][] iStateCounts = rawCombs.getCounts();
                double[] finalPs = rawCombs.getProbs();
                double[] splitPs = rawCombs.getSplitProbs();

                double sum = 0.0;
                System.out.println("State " + s + " (d,k,m)=(" + d + "," + k + "," + m + ")");
                for (int iFiss = 0; iFiss < iStateCounts.length; iFiss++) {
                    int k0 = iStateCounts[iFiss][0][1];
                    int m0 = iStateCounts[iFiss][0][2];
                    int k1 = iStateCounts[iFiss][1][1];
                    int m1 = iStateCounts[iFiss][1][2];
                    int s0 = (k0 >= 0 && m0 >= 0 && k0 + m0 <= S) ? varState[k0][m0] : -1;
                    int s1 = (k1 >= 0 && m1 >= 0 && k1 + m1 <= S) ? varState[k1][m1] : -1;
                    double p = finalPs[iFiss];
                    double pSplit = splitPs == null ? Double.NaN : splitPs[iFiss];
                    sum += p;
                    String a;
                    if (s0 >= 0 && s0 < stateCount) {
                        int k0s = kByState[s0];
                        int m0s = mByState[s0];
                        int d0 = S - k0s - m0s;
                        a = "A=" + s0 + " (" + d0 + "," + k0 + "," + m0 + ")";
                    } else {
                        a = "A=INVALID(" + s0 + ")";
                    }
                    String b;
                    if (s1 >= 0 && s1 < stateCount) {
                        int k1s = kByState[s1];
                        int m1s = mByState[s1];
                        int d1 = S - k1s - m1s;
                        b = "B=" + s1 + " (" + d1 + "," + k1 + "," + m1 + ")";
                    } else {
                        b = "B=INVALID(" + s1 + ")";
                    }
                    System.out.printf("  p=%.12g  pSplit=%.12g  %s  %s%n", p, pSplit, a, b);
                }
                System.out.printf("  sum(p)=%.15g%n", sum);
                System.out.println();
            }
        }
    }

    public static void main(String[] args) {
        if (args.length < 1 || args.length == 2 || args.length > 3) {
            System.out.println("Usage: SplitFissionCenancestorLikelihoodCore <S> [k m]");
            System.out.println("  If k and m are provided, prints only that parent state.");
            return;
        }
        int S = Integer.parseInt(args[0]);
        Integer k = null;
        Integer m = null;
        if (args.length == 3) {
            k = Integer.parseInt(args[1]);
            m = Integer.parseInt(args[2]);
        }
        printOddSplitFissionTables(S, k, m);
    }
}
