package dr.math;

//Heavily inspired in chatGPT 5

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class ThreeVHypergeometricDistribution {

    private static final int MAX_N = 200; // as before
    private static final double precission = 1e-14;

    private final int K1, K2, K3;
    private final int N;
    private final int n;

    // enumerated outcomes
    private final int[][] counts;        // [m][3]
    private final BigInteger[] numerators; // exact numerator for each outcome
    private final BigInteger denomExact;   // exact denominator = C(N,n)
    private final double[] logProbs;       // ln(prob) computed via log-factorials

    // precomputed
    private final BigInteger[][] C;
    private final double[] logFact;

    public ThreeVHypergeometricDistribution(int K1, int K2, int K3, int n) {
        this.K1 = nonNeg(K1);
        this.K2 = nonNeg(K2);
        this.K3 = nonNeg(K3);
        this.N = K1 + K2 + K3;
        this.n = n;

        if (N > MAX_N) throw new IllegalArgumentException("Total N = " + N + " exceeds MAX_N = " + MAX_N);
        if (n < 0 || n > N) throw new IllegalArgumentException("Draw size n out of range: " + n);

        this.C = precomputeBinomials(N);
        this.logFact = precomputeLogFactorials(N);

        this.denomExact = C[N][n];
        double logDenom = logBinomialDouble(N, n, logFact);

        List<int[]> triples = new ArrayList<>();
        List<BigInteger> nums = new ArrayList<>();
        List<Double> lps = new ArrayList<>();

        int x1min = Math.max(0, n - (K2 + K3));
        int x1max = Math.min(K1, n);

        for (int x1 = x1min; x1 <= x1max; x1++) {
            int remainAfter1 = n - x1;
            int x2min = Math.max(0, remainAfter1 - K3);
            int x2max = Math.min(K2, remainAfter1);
            for (int x2 = x2min; x2 <= x2max; x2++) {
                int x3 = remainAfter1 - x2;
                if (x3 < 0 || x3 > K3) continue;

                BigInteger numExact = C[K1][x1].multiply(C[K2][x2]).multiply(C[K3][x3]);

                double logNum = logBinomialDouble(K1, x1, logFact)
                        + logBinomialDouble(K2, x2, logFact)
                        + logBinomialDouble(K3, x3, logFact);
                double logP = logNum - logDenom;

                triples.add(new int[]{x1, x2, x3});
                nums.add(numExact);
                lps.add(logP);
            }
        }

        int m = triples.size();
        this.counts = new int[m][3];
        this.numerators = new BigInteger[m];
        this.logProbs = new double[m];

        for (int i = 0; i < m; i++) {
            System.arraycopy(triples.get(i), 0, this.counts[i], 0, 3);
            this.numerators[i] = nums.get(i);
            this.logProbs[i] = lps.get(i);
        }
    }

    // ACCESSORS

    /** deep copy to avoid external mutation */
    public int[][] getCounts() { return deepCopyCounts(counts); }

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

    /** exact numerators array */
    public BigInteger[] getNumerators() { return numerators.clone(); }

    /** exact denominator */
    public BigInteger getDenominator() { return denomExact; }

    // HELPERS

    private static int nonNeg(int v) {
        if (v < 0) throw new IllegalArgumentException("Count must be nonnegative: " + v);
        return v;
    }

    private static BigInteger[][] precomputeBinomials(int maxN) {
        BigInteger[][] C = new BigInteger[maxN + 1][maxN + 1];
        for (int i = 0; i <= maxN; i++) {
            for (int j = 0; j <= i; j++) {
                if (j == 0 || j == i) C[i][j] = BigInteger.ONE;
                else C[i][j] = C[i-1][j-1].add(C[i-1][j]);
            }
        }
        return C;
    }

    private static double[] precomputeLogFactorials(int maxN) {
        double[] logFact = new double[maxN + 1];
        logFact[0] = 0.0;
        for (int i = 1; i <= maxN; i++) logFact[i] = logFact[i-1] + Math.log(i);
        return logFact;
    }

    private static double logBinomialDouble(int n, int k, double[] logFact) {
        if (k < 0 || k > n) return Double.NEGATIVE_INFINITY;
        return logFact[n] - logFact[k] - logFact[n-k];
    }

    private static int[][] deepCopyCounts(int[][] src) {
        int[][] dst = new int[src.length][3];
        for (int i = 0; i < src.length; i++) System.arraycopy(src[i], 0, dst[i], 0, 3);
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

    // demo main for sanity
    public static void main(String[] args) {
        ThreeVHypergeometricDistribution d = new ThreeVHypergeometricDistribution(2, 2, 2, 6);
        int[][] counts = d.getCounts();
        double[] probs = d.getProbs();
        double[] lps = d.getLogProbs();
        BigInteger[] nums = d.getNumerators();
        BigInteger den = d.getDenominator();

        System.out.println("Den = " + den);
        double sum = 0.0;
        for (int i = 0; i < counts.length; i++) {
            System.out.printf("(%d,%d,%d)  num=%s  p=%.15g  ln p=%.8g%n",
                    counts[i][0], counts[i][1], counts[i][2],
                    nums[i].toString(), probs[i], lps[i]);
            sum += probs[i];
        }
        System.out.println("Sum p = " + sum);
    }

    @Override public boolean equals(Object o) {
        if (!(o instanceof ThreeVHypergeometricDistribution)) return false;
        ThreeVHypergeometricDistribution that = (ThreeVHypergeometricDistribution) o;
        return this.K1==that.K1 && this.K2==that.K2 && this.K3==that.K3 && this.n==that.n;
    }

    @Override public int hashCode() { return Objects.hash(K1,K2,K3,n); }
}