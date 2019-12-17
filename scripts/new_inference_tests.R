alpha <- 0
log_lambda <- 3.6
nsim <- 1000
radius <- 0.1
model <- "Geyer"
steps <- 100000
window_side <- c(0, 1)
window <- Rectangle_window(window_side, window_side)

Z <- rmultigibbs(window, matrix(alpha), exp(log_lambda), nsim = nsim, radius = radius, model = model, steps = steps)
#Z <- rppp(window, exp(log_lambda), nsim = nsim)
W <- lapply(Z, function(z) coef(gibbsm(z, window, model = model, radius = radius, print = FALSE)))
cat("Estimated values are: ", paste0(Reduce("+", W) / length(W), collapse = ", "), ".\n", sep = "")
cat("True values are: log_lambda = ", log_lambda, " and alpha = ", alpha, ".\n", sep = "")
