/* chiplot: core numerical routines (no R API, so they can be tested standalone) */
#ifndef CHIPLOT_CHICORE_H
#define CHIPLOT_CHICORE_H

void chi_biv(const double *x, const double *y, int n,
             double *lam, double *chi, int *idx, int *npts);

void chi_mult_direct(const double *X, int n, int p, const int *D, int k,
                     double *chi, double *Xrow, int *blkcount);

void chi_mult(const double *X, int n, int p, const int *D, int k,
              double *chi, double *Xrow, int *blkcount);

#endif
