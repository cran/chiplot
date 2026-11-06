/* chiplot: .Call interface and routine registration */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>
#include "chicore.h"

/* Bivariate chi-plot values: returns list(lambda, chi, index) */
SEXP C_chivals2(SEXP x_, SEXP y_)
{
    int n = LENGTH(x_);
    if (LENGTH(y_) != n) error("'x' and 'y' must have the same length");
    if (n < 2) error("need at least two observations");
    SEXP x = PROTECT(coerceVector(x_, REALSXP));
    SEXP y = PROTECT(coerceVector(y_, REALSXP));
    double *lam = (double *) R_alloc(n, sizeof(double));
    double *chi = (double *) R_alloc(n, sizeof(double));
    int *idx = (int *) R_alloc(n, sizeof(int));
    int npts = 0;

    chi_biv(REAL(x), REAL(y), n, lam, chi, idx, &npts);

    SEXP rlam = PROTECT(allocVector(REALSXP, npts));
    SEXP rchi = PROTECT(allocVector(REALSXP, npts));
    SEXP ridx = PROTECT(allocVector(INTSXP, npts));
    for (int i = 0; i < npts; i++) {
        REAL(rlam)[i] = lam[i];
        REAL(rchi)[i] = chi[i];
        INTEGER(ridx)[i] = idx[i];
    }
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_VECTOR_ELT(out, 0, rlam); SET_STRING_ELT(nms, 0, mkChar("lambda"));
    SET_VECTOR_ELT(out, 1, rchi); SET_STRING_ELT(nms, 1, mkChar("chi"));
    SET_VECTOR_ELT(out, 2, ridx); SET_STRING_ELT(nms, 2, mkChar("index"));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(7);
    return out;
}

/* Multivariate chi-values: X is an n x p numeric matrix, D integer block sizes */
SEXP C_chimult(SEXP X_, SEXP D_)
{
    SEXP dim = getAttrib(X_, R_DimSymbol);
    if (!isMatrix(X_) || LENGTH(dim) != 2) error("'X' must be a matrix");
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    SEXP X = PROTECT(coerceVector(X_, REALSXP));
    SEXP D = PROTECT(coerceVector(D_, INTSXP));
    int k = LENGTH(D), tot = 0;
    for (int l = 0; l < k; l++) {
        if (INTEGER(D)[l] < 1) error("all block sizes in 'D' must be positive");
        tot += INTEGER(D)[l];
    }
    if (tot != p) error("block sizes in 'D' must sum to the number of columns");

    SEXP chi = PROTECT(allocVector(REALSXP, n));
    double *Xrow = (double *) R_alloc((size_t) n * p, sizeof(double));
    int *blk = (int *) R_alloc(k, sizeof(int));
    chi_mult(REAL(X), n, p, INTEGER(D), k, REAL(chi), Xrow, blk);
    UNPROTECT(3);
    return chi;
}

static const R_CallMethodDef CallEntries[] = {
    {"C_chivals2", (DL_FUNC) &C_chivals2, 2},
    {"C_chimult",  (DL_FUNC) &C_chimult,  2},
    {NULL, NULL, 0}
};

void R_init_chiplot(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
