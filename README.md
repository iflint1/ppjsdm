# Saturated pairwise interaction Gibbs point process model fitting and simulation

`ppjsdm` is an R package for working with different types of saturated pairwise interaction Gibbs point processes.
It is intended to be used in the context of joint species distribution models, i.e. datasets containing the locations of individuals along with marks representing their species.
The package relies on parts of `spatstat`, in particular their `im` objects which in our framework can be used as windows or covariates.

The technical details of the point process model are presented in this paper: https://doi.org/10.1111/rssc.12596. 
The package was presented at ISEC 2020, see the `isec2020_presentation_handout` vignette for the slides that show the package in action.


## Installation Requirements
For package installation, gfortan and a c++ complier are necessary. Therefore, the R environment should be set up as following: 

For Windows users, the `Rtools` toolchain is required, see [the following tutorial](https://cran.r-project.org/bin/windows/Rtools/) for help installing `Rtools`. Check RTools version is compatible R version (i.e. RTools 4.5 = R 4.5)

Mac users require xcode and a c++ complier to install the package. See the following guides: for [xcode and gfortran](https://mac.r-project.org/tools/), or for [gfortran downloading help](https://cran.r-project.org/bin/macosx/tools/). For recent versions of Mac and R complier 'gfortran-12.2-universal.pkg' required. For other versions of Mac and R, check version compatibility in the above guides. 


## Tutorials 
An in-depth tutorial on model application including parametrisation, fitting, prediction and troubleshooting is available at <https://github.com/shenali-fernando/ppjsdm_tutorial/>. More tutorials can be found in the repo. 

The accompanying package `ppjsdm_on_datasets` demonstrates the package on a few datasets from plant ecology.


**Note**: The `ppjsdm` package is still in active development, so please do not expect any kind of stability of the user interface.



<!-- badges: start -->
  [![R build status](https://github.com/iflint1/ppjsdm/workflows/R-CMD-check/badge.svg)](https://github.com/iflint1/ppjsdm/actions)
<!-- badges: end -->
