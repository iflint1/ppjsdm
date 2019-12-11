alpha <- -1
log_lambda <- 3.6
nsim <- 1000
radius <- 0.1
model <- "Geyer"

Z <- rmultigibbs(Rectangle_window(), matrix(alpha), exp(log_lambda), nsim = nsim, radius = radius, model = model)
W <- lapply(Z, function(z) coef(gibbsm(z, Rectangle_window(), model = model, radius = radius, print = FALSE)))
cat("Estimated values are: ", paste0(Reduce("+", W) / length(W), collapse = ", "), ".\n", sep = "")
cat("True values are: log_lambda = ", log_lambda, " and alpha = ", alpha, ".\n", sep = "")
