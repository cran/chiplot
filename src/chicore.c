/*
 * chiplot: core numerical routines for chi-plots.
 *
 * chi_biv  -- bivariate chi-plot of Fisher & Switzer (1985, 2001).
 *             A direct C translation of the Fortran subroutine CHIPLOT (chisub.f).
 * chi_mult -- multivariate chi-values of Fisher & Lee (2026), using generalized
 *             orthants. Gives the same values as the Rcpp function chivalsC
 *             (Hval.cpp), but evaluates the full orthant and all block orthants
 *             in a single pass over the data.
 */
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "chicore.h"

#if defined(__GNUC__) || defined(__clang__)
#define POPCNT64(w) __builtin_popcountll(w)
#else
static int POPCNT64(uint64_t w)
{
    w = w - ((w >> 1) & 0x5555555555555555ULL);
    w = (w & 0x3333333333333333ULL) + ((w >> 2) & 0x3333333333333333ULL);
    w = (w + (w >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (int) ((w * 0x0101010101010101ULL) >> 56);
}
#endif

/*
 * Bivariate chi-plot.
 *   x, y : data of length n (no missing values)
 *   lam  : output lambda values (length >= n)
 *   chi  : output chi values    (length >= n)
 *   idx  : output 1-based indices of the points retained (length >= n)
 *   npts : number of points retained
 * For each point, F, G and H are the marginal and joint empirical distribution
 * functions evaluated at that point, computed from the other n - 1 points.
 * Points with |lambda| >= 4 (1/(n-1) - 1/2)^2, or with F, G in {0, 1}, are dropped.
 */
void chi_biv(const double *x, const double *y, int n,
             double *lam, double *chi, int *idx, int *npts)
{
    double nm1 = (double)(n - 1);
    double xlim = 2.0 * (1.0 / nm1 - 0.5);
    xlim = xlim * xlim;
    int ncalc = 0;

    for (int i = 0; i < n; i++) {
        double xi = x[i], yi = y[i];
        int mx = 0, my = 0, mh = 0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            int a = (x[j] <= xi), b = (y[j] <= yi);
            mx += a;
            my += b;
            mh += (a & b);
        }
        double F = (double)mx / nm1, G = (double)my / nm1, H = (double)mh / nm1;
        double f5 = F - 0.5, g5 = G - 0.5;
        double sgn = (f5 * g5 < 0.0) ? -1.0 : 1.0;
        double m2 = (f5 * f5 > g5 * g5) ? f5 * f5 : g5 * g5;
        double tit = 4.0 * sgn * m2;
        if (fabs(tit) >= xlim) continue;
        double uv = F * G * (1.0 - F) * (1.0 - G);
        if (uv == 0.0) continue;
        lam[ncalc] = tit;
        chi[ncalc] = (H - F * G) / sqrt(uv);
        idx[ncalc] = i + 1;
        ncalc++;
    }
    *npts = ncalc;
}

/*
 * Multivariate chi-values.
 *   X        : n x p data matrix, column-major (columns assumed median-centred)
 *   D        : block sizes D[0..k-1], summing to p; blocks are consecutive columns
 *   chi      : output, length n
 *   Xrow     : workspace, length n * p (row-major copy of X)
 *   blkcount : workspace, length k
 *
 * For data point i, the generalized orthant is
 *   { v : v_c <= x_ic  if x_ic > 0;  v_c > x_ic  if x_ic <= 0 },  c = 1..p.
 * H_i is the proportion of the n data points in the orthant and H_li the
 * proportion in the orthant restricted to the columns of block l. Proportions
 * below 1e-8 are set to 1/n and those above 1 - 1e-8 to 1 - 1/n (as in Hval.cpp).
 * Then
 *   chi_i = sqrt(n) (H_i - prod_l H_li) / sqrt( P (1 - P (1 + sum_l (1 - H_li)/H_li)) ),
 * with P = prod_l H_li.
 */
/* Turn orthant counts into chi-values (shared by both algorithms) */
static double chi_from_counts(int full, const int *blkcount, int k, int n)
{
    double dn = (double) n;
    double hval = (double) full / dn;
    if (hval < 1e-8) hval = 1.0 / dn;
    if (hval > 1.0 - 1e-8) hval = 1.0 - 1.0 / dn;
    double sdSum = 1.0, prod = 1.0;
    for (int l = 0; l < k; l++) {
        double h = (double) blkcount[l] / dn;
        if (h < 1e-8) h = 1.0 / dn;
        if (h > 1.0 - 1e-8) h = 1.0 - 1.0 / dn;
        sdSum += (1.0 - h) / h;
        prod *= h;
    }
    return sqrt(dn) * (hval - prod) / sqrt(prod * (1.0 - prod * sdSum));
}

/* Direct algorithm: O(n^2 p) comparisons, O(n p) memory */
void chi_mult_direct(const double *X, int n, int p, const int *D, int k,
                     double *chi, double *Xrow, int *blkcount)
{
    for (int i = 0; i < n; i++)
        for (int c = 0; c < p; c++)
            Xrow[(long) i * p + c] = X[(long) c * n + i];

    for (int i = 0; i < n; i++) {
        const double *xi = Xrow + (long) i * p;
        int full = 0;
        for (int l = 0; l < k; l++) blkcount[l] = 0;

        for (int j = 0; j < n; j++) {
            const double *xj = Xrow + (long) j * p;
            int allin = 1, c = 0;
            for (int l = 0; l < k; l++) {
                int inb = 1, cend = c + D[l];
                for (; c < cend; c++) {
                    int in = (xi[c] > 0.0) ? (xj[c] <= xi[c]) : (xj[c] > xi[c]);
                    if (!in) { inb = 0; c = cend; break; }
                }
                c = cend;
                if (inb) blkcount[l]++; else allin = 0;
            }
            full += allin;
        }

        chi[i] = chi_from_counts(full, blkcount, k, n);
    }
}

/* ------------------------------------------------------------------------
 * Bit-set algorithm: O(n log n p + n^2 p / 64) operations, O(n^2 / 32) bytes.
 * For column c, the set { j : x_jc <= v } is a prefix of the rows sorted by
 * column c. Running through the sorted rows while building that prefix as a
 * bit set, each point's orthant is updated with one AND (or AND NOT) per
 * 64 rows. Counts, and hence chi-values, are identical to chi_mult_direct.
 * ------------------------------------------------------------------------ */
typedef struct { double v; int i; } vi_pair;

static int cmp_vi(const void *a, const void *b)
{
    const vi_pair *x = (const vi_pair *) a, *y = (const vi_pair *) b;
    if (x->v < y->v) return -1;
    if (x->v > y->v) return 1;
    return (x->i > y->i) - (x->i < y->i);
}

/* work sizes: bits must hold 2 * n * W + W words, W = (n + 63) / 64 */
static int chi_mult_bits(const double *X, int n, const int *D, int k,
                  double *chi, int *blkcount_all, uint64_t *bits, vi_pair *pairs)
{
    int W = (n + 63) / 64;
    uint64_t *acc = bits;                       /* n x W, block accumulators */
    uint64_t *full = bits + (size_t) n * W;     /* n x W, full orthant       */
    uint64_t *cur = full + (size_t) n * W;      /* W, current prefix         */
    uint64_t lastmask = (n % 64 == 0) ? ~(uint64_t) 0 : (((uint64_t) 1 << (n % 64)) - 1);

    for (size_t t = 0; t < (size_t) n * W; t++) full[t] = ~(uint64_t) 0;
    for (int i = 0; i < n; i++) full[(size_t) i * W + W - 1] &= lastmask;

    int c = 0;
    for (int l = 0; l < k; l++) {
        for (size_t t = 0; t < (size_t) n * W; t++) acc[t] = ~(uint64_t) 0;
        for (int i = 0; i < n; i++) acc[(size_t) i * W + W - 1] &= lastmask;

        for (int cend = c + D[l]; c < cend; c++) {
            const double *col = X + (size_t) c * n;
            for (int i = 0; i < n; i++) { pairs[i].v = col[i]; pairs[i].i = i; }
            qsort(pairs, (size_t) n, sizeof(vi_pair), cmp_vi);
            for (int w = 0; w < W; w++) cur[w] = 0;
            int start = 0;
            for (int pos = 0; pos < n; pos++) {
                int j = pairs[pos].i;
                cur[j >> 6] |= (uint64_t) 1 << (j & 63);
                if (pos == n - 1 || pairs[pos + 1].v != pairs[pos].v) {
                    /* tie group start..pos: all share the set { x <= v } = cur */
                    for (int q = start; q <= pos; q++) {
                        int i = pairs[q].i;
                        uint64_t *a = acc + (size_t) i * W;
                        if (col[i] > 0.0) {
                            for (int w = 0; w < W; w++) a[w] &= cur[w];
                        } else {
                            for (int w = 0; w < W; w++) a[w] &= ~cur[w];
                        }
                    }
                    start = pos + 1;
                }
            }
        }
        for (int i = 0; i < n; i++) {
            uint64_t *a = acc + (size_t) i * W, *f = full + (size_t) i * W;
            int cnt = 0;
            for (int w = 0; w < W; w++) { cnt += POPCNT64(a[w]); f[w] &= a[w]; }
            blkcount_all[(size_t) l * n + i] = cnt;
        }
    }

    int *bc = (int *) malloc((size_t) k * sizeof(int));
    if (!bc) return 1;
    for (int i = 0; i < n; i++) {
        uint64_t *f = full + (size_t) i * W;
        int cnt = 0;
        for (int w = 0; w < W; w++) cnt += POPCNT64(f[w]);
        for (int l = 0; l < k; l++) bc[l] = blkcount_all[(size_t) l * n + i];
        chi[i] = chi_from_counts(cnt, bc, k, n);
    }
    free(bc);
    return 0;
}

/* Dispatcher: bit-set algorithm unless its workspace would exceed ~1 GB */
void chi_mult(const double *X, int n, int p, const int *D, int k,
              double *chi, double *Xrow, int *blkcount)
{
    size_t W = (size_t) (n + 63) / 64;
    size_t nbits = 2 * (size_t) n * W + W;
    if (nbits * sizeof(uint64_t) <= ((size_t) 1 << 30)) {
        uint64_t *bits = (uint64_t *) malloc(nbits * sizeof(uint64_t));
        vi_pair *pairs = (vi_pair *) malloc((size_t) n * sizeof(vi_pair));
        int *bca = (int *) malloc((size_t) n * k * sizeof(int));
        int fail = (!bits || !pairs || !bca);
        if (!fail) fail = chi_mult_bits(X, n, D, k, chi, bca, bits, pairs);
        free(bits); free(pairs); free(bca);
        if (!fail) return;
    }
    chi_mult_direct(X, n, p, D, k, chi, Xrow, blkcount);
}
