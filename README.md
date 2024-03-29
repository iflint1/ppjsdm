# Saturated pairwise interaction Gibbs point process model fitting and simulation

`ppjsdm` is an R package for working with different types of saturated pairwise interaction Gibbs point processes.
It is intended to be used in the context of joint species distribution models, i.e. datasets containing the locations of individuals along with marks representing their species.
The package relies on parts of `spatstat`, in particular their `im` objects which in our framework can be used as windows or covariates.

The package was presented at ISEC 2020, see the `isec2020_presentation_handout` vignette for the slides that show the package in action.

The accompanying package `ppjsdm_on_datasets` demonstrates the package on a few datasets from plant ecology.

The `Rtools` toolchain is required to build the package on Windows, see [the following tutorial](https://cran.r-project.org/bin/windows/Rtools/) for help installing `Rtools`.

**Note**: The `ppjsdm` package is still in active development, so please do not expect any kind of stability of the user interface.

**Note**: The technical details of the point process model will be presented in a forthcoming paper, which will be linked to here.

<!-- badges: start -->
  [![R build status](https://github.com/iflint1/ppjsdm/workflows/R-CMD-check/badge.svg)](https://github.com/iflint1/ppjsdm/actions)
<!-- badges: end -->
