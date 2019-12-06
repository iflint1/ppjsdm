p <- 2

repulsion <- 10
alpha <- diag(p) * -repulsion + repulsion
print(alpha)

lambda <- c(100, 100)
print(lambda)

steps <- 400
sim_per_step <- 10
pause <- 0.05

chain <- rgibbs(alpha, lambda, nsim = 1 + (steps - 1) * sim_per_step, return_full_chain = TRUE)

for(i in seq_len(steps)) {
  plot(chain[[1 + (i - 1) * sim_per_step]], Rectangle_window())
  Sys.sleep(pause)
}
