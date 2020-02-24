library(ppjsdm)
library(spatstat)

N <- 10000
confidence <- 0.95
window <- Rectangle_window(c(0, 2), c(-1, 1))
spatstat_window <- owin(c(0, 2), c(-1, 1))

lambda <- 30
radius <- 0.1

ppjsdm_configurations <- rgibbs(window = window, alpha = matrix(-Inf), lambda = lambda, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rHardcore(lambda, W = spatstat_window, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
       "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))


lambda <- 35
radius <- 0.05
gamma <- 0.6

ppjsdm_configurations <- rgibbs(window = window, alpha = matrix(log(gamma)), lambda = lambda, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rStrauss(lambda, W = spatstat_window, gamma = gamma, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
               "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))
