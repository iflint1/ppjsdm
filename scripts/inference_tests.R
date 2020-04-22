library(ppjsdm)

alpha <- matrix(-0.1)
gamma <- matrix(0)
lambda <- 30
nsim <- 1000
short_range <- matrix(0.1)
medium_range <- matrix(0)
long_range <- matrix(0)
model <- "square_exponential"
medium_range_model <- "square_exponential"
steps <- 0
window_side <- c(0, 1)
window <- Rectangle_window(window_side, window_side)
use_rgibbs <- TRUE
use_joint_pdf <- TRUE

if(use_rgibbs) {
  Z <- ppjsdm::rgibbs(window = window,
                      alpha = alpha,
                      gamma = gamma,
                      lambda = lambda,
                      nsim = nsim,
                      short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      model = model,
                      medium_range_model = medium_range_model,
                      steps = steps,
                      drop = FALSE)
} else {
  Z <- ppjsdm::rppp(window = window,
                    lambda = lambda,
                    nsim = nsim,
                    drop = FALSE)
}

if(use_joint_pdf) {
  G <-gibbsm(Z,
             window = window,
             model = model,
             medium_range_model = medium_range_model,
             short_range = short_range,
             medium_range = medium_range,
             long_range = long_range,
             print = FALSE,
             use_glmnet = FALSE)
  estimate <- G$coefficients
} else {
  G <-lapply(Z, function(z) gibbsm(z,
                                   window = window,
                                   model = model,
                                   medium_range_model = medium_range_model,
                                   short_range = short_range,
                                   medium_range = medium_range,
                                   long_range = long_range,
                                   print = FALSE,
                                   use_glmnet = FALSE))
  W <- lapply(G, function(g) g$coefficients)
  if(length(W) == 1) {
    estimate <- unlist(W)
  } else {
    estimate <- Reduce(function(a, b) unlist(a) + unlist(b), W) / length(W)
  }
}
cat("Estimated values are: ", paste0(estimate, collapse = ", "), ".\n", sep = "")
cat("True values are: lambda = ", lambda, " and alpha = ", alpha, ".\n", sep = "")
