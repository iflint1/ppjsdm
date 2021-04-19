library(microbenchmark)
library(ppjsdm)
library(spatstat)
remove(list = ls())

window_spatstat <- owin()

x_covariate <- function(x, y) x
y_covariate <- function(x, y) y

dummy_factor <- 1e6
model <- "Geyer"
medium_range_model <- "Geyer"

set.seed(42)

configuration_spatstat <- rpoispp(lambda = function(x, y) exp(5 + x_covariate(x, y) + y_covariate(x, y)), win = window_spatstat)
configuration <- as.Configuration(configuration_spatstat)

plot(configuration_spatstat)

nd <- 100

b <- microbenchmark(
  "ppjsdm::gibbsm with short-range and medium-range" = ppjsdm::gibbsm(configuration,
                                                                      covariates = list(x_covariate = x_covariate, y_covariate = y_covariate),
                                                                      use_glmnet = FALSE,
                                                                      short_range = 0.01,
                                                                      medium_range = 0.02,
                                                                      long_range = 0.04,
                                                                      dummy_factor = dummy_factor,
                                                                      max_dummy = nd * nd,
                                                                      model = model,
                                                                      medium_range_model = medium_range_model),
  "ppjsdm::gibbsm with short-range" = ppjsdm::gibbsm(configuration,
                                                     covariates = list(x_covariate = x_covariate, y_covariate = y_covariate),
                                                     use_glmnet = FALSE,
                                                     short_range = 0.01,
                                                     dummy_factor = dummy_factor,
                                                     max_dummy = nd * nd,
                                                     model = model,
                                                     medium_range_model = medium_range_model),
  "ppjsdm::gibbsm without ranges" = ppjsdm::gibbsm(configuration,
                                                   covariates = list(x_covariate = x_covariate, y_covariate = y_covariate),
                                                   use_glmnet = FALSE,
                                                   short_range = 0,
                                                   dummy_factor = dummy_factor,
                                                   max_dummy = nd * nd,
                                                   model = model,
                                                   medium_range_model = medium_range_model),
  # `logi` is the fastest method according to docs, but `mpl` appears to be faster here.
  "spatstat::ppm Poisson mpl" = ppm(configuration_spatstat ~ 1 + x_covariate + y_covariate, method = "mpl", nd = c(nd, nd)),
  "spatstat::ppm Poisson logi" = ppm(configuration_spatstat ~ 1 + x_covariate + y_covariate, method = "logi", nd = c(nd, nd)),
  "spatstat::ppm Strauss mpl" = ppm(configuration_spatstat ~ 1 + x_covariate + y_covariate, method = "mpl", nd = c(nd, nd), correction = "none", interaction = Strauss(0.01)),
  "spatstat::ppm Strauss logi" = ppm(configuration_spatstat ~ 1 + x_covariate + y_covariate, method = "logi", nd = c(nd, nd), correction = "none", interaction = Strauss(0.01)),
  times = 10L
)

print(b)
