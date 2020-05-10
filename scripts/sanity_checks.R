remove(list = ls())
library(ppjsdm)
library(spatstat)

N <- 10000
confidence <- 0.95
window <- Rectangle_window(c(0, 2), c(-1, 1))
spatstat_window <- owin(window$x_range, window$y_range)

beta0 <- log(30)
radius <- 0.1

ppjsdm_configurations <- rgibbs(window = window, alpha = matrix(-Inf), beta0 = beta0, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rHardcore(exp(beta0), W = spatstat_window, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
       "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))


beta0 <- log(35)
radius <- 0.05
gamma <- 0.6

ppjsdm_configurations <- rgibbs(window = window, alpha = matrix(0.5 * log(gamma)), beta0 = beta0, gamma = 0, saturation = Inf, model = "Geyer", short_range = radius, nsim = N)
ppjsdm_number_points <- sapply(ppjsdm_configurations, function(a) length(a$x))
ppjsdm_result <- t.test(ppjsdm_number_points, conf.interval = confidence)

spatstat_configurations <- rStrauss(exp(beta0), W = spatstat_window, gamma = gamma, R = radius, expand = FALSE, nsim = N)
spatstat_number_points <- sapply(spatstat_configurations, function(a) length(a$x))
spatstat_result <- t.test(spatstat_number_points, conf.interval = confidence)

message(paste0("ppjsdm confidence interval: [", paste0(ppjsdm_result$conf.int, collapse = ", "),
               "], spatstat confidence interval: [", paste0(spatstat_result$conf.int, collapse = ", "), "]."))

window <- Rectangle_window(c(0, 1), c(0, 1))
alpha <- 0.
beta0 <- log(50)
gamma <- -0.4
model <- "Geyer"
medium_range_model <- "Geyer"
short_range <- 0.05
medium_range <- 0.1
long_range <- 0.15
saturation <- 2


N <- 10000
cftp_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N)
cftp_number_points <- sapply(cftp_configurations, function(a) length(a$x))
cftp_result <- t.test(cftp_number_points, conf.interval = confidence)

N <- 1000
steps <- 100000
mh_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N, steps = steps)
mh_number_points <- sapply(mh_configurations, function(a) length(a$x))
mh_result <- t.test(mh_number_points, conf.interval = confidence)

message(paste0("cftp confidence interval: [", paste0(cftp_result$conf.int, collapse = ", "),
               "], mh confidence interval: [", paste0(mh_result$conf.int, collapse = ", "), "]."))

window <- Rectangle_window(c(0, 1), c(0, 1))
alpha <- cbind(c(-1.0, 0.), c(0., 0.))
beta0 <- c(log(50), -Inf)
gamma <- cbind(c(0., 0.), c(0., 0.))
model <- "Geyer"
medium_range_model <- "Geyer"
short_range <- matrix(0.05, 2, 2)
medium_range <- matrix(0.1, 2, 2)
long_range <- matrix(0.15, 2, 2)
saturation <- 2


N <- 10000
cftp_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N)
cftp_number_points <- sapply(cftp_configurations, function(a) length(a$x))
cftp_result <- t.test(cftp_number_points, conf.interval = confidence)

N <- 1000
steps <- 100000
mh_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N, steps = steps)
mh_number_points <- sapply(mh_configurations, function(a) length(a$x))
mh_result <- t.test(mh_number_points, conf.interval = confidence)

message(paste0("cftp confidence interval: [", paste0(cftp_result$conf.int, collapse = ", "),
               "], mh confidence interval: [", paste0(mh_result$conf.int, collapse = ", "), "]."))

window <- Rectangle_window(c(0, 1), c(0, 1))
alpha <- cbind(c(-1.0, -0.3), c(-0.3, -0.4))
beta0 <- c(log(30), log(30))
gamma <- cbind(c(0, 0), c(0, 0))
model <- "square_bump"
medium_range_model <- "square_exponential"
short_range <- matrix(0.05, 2, 2)
medium_range <- matrix(0.1, 2, 2)
long_range <- matrix(0.15, 2, 2)
saturation <- 2


N <- 10000
cftp_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N)
cftp_number_points <- sapply(cftp_configurations, function(a) length(a$x))
cftp_result <- t.test(cftp_number_points, conf.interval = confidence)

N <- 1000
steps <- 100000
mh_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N, steps = steps)
mh_number_points <- sapply(mh_configurations, function(a) length(a$x))
mh_result <- t.test(mh_number_points, conf.interval = confidence)

message(paste0("cftp confidence interval: [", paste0(cftp_result$conf.int, collapse = ", "),
               "], mh confidence interval: [", paste0(mh_result$conf.int, collapse = ", "), "]."))

window <- Rectangle_window(c(0, 1), c(0, 1))
alpha <- 0.4
beta0 <- log(10)
gamma <- 0
model <- "Geyer"
medium_range_model <- "Geyer"
short_range <- 0.05
medium_range <- 0.1
long_range <- 0.15
saturation <- 1


N <- 10000
cftp_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N)
cftp_number_points <- sapply(cftp_configurations, function(a) length(a$x))
cftp_result <- t.test(cftp_number_points, conf.interval = confidence)

N <- 1000
steps <- 100000
mh_configurations <- rgibbs(window = window, saturation = saturation, alpha = alpha, beta0 = beta0, gamma = gamma, model = model, medium_range_model = medium_range_model, short_range = short_range, medium_range = medium_range, long_range = long_range, nsim = N, steps = steps)
mh_number_points <- sapply(mh_configurations, function(a) length(a$x))
mh_result <- t.test(mh_number_points, conf.interval = confidence)

message(paste0("cftp confidence interval: [", paste0(cftp_result$conf.int, collapse = ", "),
               "], mh confidence interval: [", paste0(mh_result$conf.int, collapse = ", "), "]."))

