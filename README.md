# Multivariate Gibbs point process model fitting and simulation

`ppjsdm` is an R package for working with different types of multivariate Gibbs point processes.
It is intended to be used in the context of joint species distribution models, i.e. datasets containing the locations of individuals along with marks representing their species.
The package is based on `spatstat`, mainly for their `im` objects which in our framework can be used as windows or covariates.
