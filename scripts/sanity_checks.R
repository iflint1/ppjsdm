library(ppjsdm)
library(spatstat)

N <- 100000
confidence <- 0.95

lambda <- 40
radius <- 0.1

ppjsdm_configurations <- rgibbs(alpha = matrix(-Inf), lambda = lambda, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rHardcore(lambda, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
       "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))


lambda <- 55
radius <- 0.05
gamma <- 0.6

ppjsdm_configurations <- rgibbs(alpha = matrix(log(gamma)), lambda = lambda, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rStrauss(lambda, gamma = gamma, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
               "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))
