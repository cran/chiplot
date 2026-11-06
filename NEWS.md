# chiplot 0.1.0

* First version. Functions `chiplot()`, `chivals()`, `chipairs()`, `chimult()`,
  `chipart()`, `chidepth()` and `chi.partitions()`, with `plot()` and `print()`
  methods.
* Compiled code is in C: the bivariate chi-values (from the Fortran subroutine
  in `chisub.f`) and the multivariate chi-values (from `Hval.cpp`, using a
  bit-set algorithm that gives identical values much faster).
* Data depths are computed in base R; `"RP"`, `"TD"` and `"MhD"` reproduce
  `fda.usc::mdepth.RP`, `mdepth.TD` and `mdepth.MhD`. With the same seed,
  `chimult()` reproduces the original `chimult.lam()` exactly.
