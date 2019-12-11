alpha <- -0.6
log_lambda <- 0.36
nsim <- 1000
radius <- 1
model <- "Geyer"

Z <- rmultigibbs(Rectangle_window(c(0, 10), c(0, 10)), matrix(alpha), exp(log_lambda), nsim = nsim, radius = radius, model = model)
W <- lapply(Z, function(z) coef(gibbsm(z, Rectangle_window(c(0, 10), c(0, 10)), model = model, radius = radius, print = FALSE)))
cat("Estimated values are: ", paste0(Reduce("+", W) / length(W), collapse = ", "), ".\n", sep = "")
cat("True values are: log_lambda = ", log_lambda, " and alpha = ", alpha, ".\n", sep = "")
