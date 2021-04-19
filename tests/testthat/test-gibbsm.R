library(ppjsdm)

context("gibbsm")

# test_that("gibbsm returns correct parameters for Geyer potentials potentials", {
#   short_range <- matrix(0.01, 2, 2)
#   medium_range <- matrix(0.02, 2, 2)
#   long_range <- matrix(0.04, 2, 2)
#   covariates <- list(x = function(x, y) x)
#   max_dummy <- 1e3
#   dummy_factor <- 1e6
#   saturation <- 2
#
#   set.seed(1)
#   z <- ppjsdm::rppp(lambda = c(1e2, 1e2))
#
#   fit <- gibbsm(z, short_range = short_range,
#                 medium_range = medium_range,
#                 long_range = long_range,
#                 use_glmnet = FALSE,
#                 max_dummy = max_dummy,
#                 covariates = covariates,
#                 dummy_factor = dummy_factor,
#                 saturation = saturation,
#                 model = "Geyer",
#                 medium_range_model = "Geyer")
#
#   target_alpha <- cbind(c(-6.59959410750887087715454981662333011627197265625000, -0.28686623020021911889543275719915982335805892944336),
#                         c(-0.28686623020021911889543275719915982335805892944336, 0.30721670699005437787931782622763421386480331420898))
#   target_gamma <- cbind(c(-0.05388042463908571805264458021156315226107835769653, 0.05478424262149115403497390275333600584417581558228),
#                         c(0.05478424262149115403497390275333600584417581558228, -0.01078111309813421098136032583170162979513406753540))
#   target_beta0 <- c(4.66835323901555110381877966574393212795257568359375, 4.73670372000843808990566685679368674755096435546875)
#
#   expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
#   expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
#   expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
# })
#
# test_that("gibbsm returns correct parameters for Exponential/Geyer potentials", {
#   short_range <- matrix(0.01, 2, 2)
#   medium_range <- matrix(0.02, 2, 2)
#   long_range <- matrix(0.04, 2, 2)
#   covariates <- list(x = function(x, y) x)
#   max_dummy <- 1e3
#   dummy_factor <- 1e6
#   saturation <- 2
#
#   set.seed(1)
#   z <- ppjsdm::rppp(lambda = c(1e2, 1e2))
#
#   fit <- gibbsm(z, short_range = short_range,
#                 medium_range = medium_range,
#                 long_range = long_range,
#                 use_glmnet = FALSE,
#                 max_dummy = max_dummy,
#                 covariates = covariates,
#                 dummy_factor = dummy_factor,
#                 saturation = saturation,
#                 model = "exponential",
#                 medium_range_model = "Geyer")
#
#   target_alpha <- cbind(c(0.36336369048743522025546326403855346143245697021484, -0.19446700039060069165053334927506512030959129333496),
#                         c(-0.19446700039060069165053334927506512030959129333496, 0.27457522734498557293036924420448485761880874633789))
#   target_gamma <- cbind(c(-0.08297137548644326066060727953299647197127342224121, 0.08329905314773396185490383913929690606892108917236),
#                         c(0.08329905314773396185490383913929690606892108917236, -0.04220028235639199098994112091531860642135143280029))
#   target_beta0 <- c(4.57327077246616475747487129410728812217712402343750, 4.72239072452669272905723119038157165050506591796875)
#
#   expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
#   expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
#   expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
# })
#
# test_that("gibbsm returns correct parameters for Exponential/Half exponential potentials", {
#   short_range <- matrix(0.01, 2, 2)
#   medium_range <- matrix(0.02, 2, 2)
#   long_range <- matrix(0.04, 2, 2)
#   covariates <- list(x = function(x, y) x)
#   max_dummy <- 1e3
#   dummy_factor <- 1e6
#   saturation <- 2
#
#   set.seed(1)
#   z <- ppjsdm::rppp(lambda = c(1e2, 1e2))
#
#   # Exponential/Half exponential potentials
#   fit <- gibbsm(z, short_range = short_range,
#                 medium_range = medium_range,
#                 long_range = long_range,
#                 use_glmnet = FALSE,
#                 max_dummy = max_dummy,
#                 covariates = covariates,
#                 dummy_factor = dummy_factor,
#                 saturation = saturation,
#                 model = "exponential",
#                 medium_range_model = "half_exponential")
#
#   target_alpha <- cbind(c(0.46316221468203372380045834688644390553236007690430, -0.12841035124783206633480858727125450968742370605469),
#                         c(-0.12841035124783206633480858727125450968742370605469, 0.26250071426500326543518326616322156041860580444336))
#   target_gamma <- cbind(c(-0.35514469623914096674255347352300304919481277465820, 0.04449345621518227544832768671767553314566612243652),
#                         c(0.04449345621518227544832768671767553314566612243652, -0.02976090796946117764121275683919520815834403038025))
#   target_beta0 <- c(4.83942274113594095297230524010956287384033203125000, 4.72348555684051074621265797759406268596649169921875)
#
#   expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
#   expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
#   expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
# })
#
# test_that("gibbsm returns correct parameters for Exponential/Square exponential potentials", {
#   short_range <- matrix(0.01, 2, 2)
#   medium_range <- matrix(0.02, 2, 2)
#   long_range <- matrix(0.04, 2, 2)
#   covariates <- list(x = function(x, y) x)
#   max_dummy <- 1e3
#   dummy_factor <- 1e6
#   saturation <- 2
#
#   set.seed(1)
#   z <- ppjsdm::rppp(lambda = c(1e2, 1e2))
#
#   fit <- gibbsm(z, short_range = short_range,
#                 medium_range = medium_range,
#                 long_range = long_range,
#                 use_glmnet = FALSE,
#                 max_dummy = max_dummy,
#                 covariates = covariates,
#                 dummy_factor = dummy_factor,
#                 saturation = saturation,
#                 model = "exponential",
#                 medium_range_model = "square_exponential")
#
#   target_alpha <- cbind(c(0.50757549646651989316836761645390652120113372802734, -0.28208707564593760164228797293617390096187591552734),
#                         c(-0.28208707564593760164228797293617390096187591552734, 0.33349752459166975837590030096180271357297897338867))
#   target_gamma <- cbind(c(-0.20753037003541913163573440215259324759244918823242, 0.14356938938182947640420650259329704567790031433105),
#                         c(0.14356938938182947640420650259329704567790031433105, -0.09404490341622943894162034439432318322360515594482))
#   target_beta0 <- c(4.59347691764322263452413608320057392120361328125000, 4.71948357859354583609956534928642213344573974609375)
#
#   expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
#   expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
#   expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
# })

test_that("gibbsm works with infinite saturation for Geyer potentials potentials", {
  ntypes <- 4
  short_range <- matrix(0.01, ntypes, ntypes)
  medium_range <- matrix(0.02, ntypes, ntypes)
  long_range <- matrix(0.04, ntypes, ntypes)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 1e6

  set.seed(1)
  z <- ppjsdm::rppp(lambda = rep(1e2, ntypes))

  set.seed(1)
  fit_large <- gibbsm(z, short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      use_glmnet = FALSE,
                      max_dummy = max_dummy,
                      dummy_factor = dummy_factor,
                      saturation = 1e6,
                      model = "Geyer",
                      medium_range_model = "Geyer")

  set.seed(1)
  fit_infinite <- gibbsm(z, short_range = short_range,
                         medium_range = medium_range,
                         long_range = long_range,
                         use_glmnet = FALSE,
                         max_dummy = max_dummy,
                         dummy_factor = dummy_factor,
                         saturation = Inf,
                         model = "Geyer",
                         medium_range_model = "Geyer")

  expect_equal(as.vector(fit_large$coefficients$alpha), as.vector(fit_infinite$coefficients$alpha))
  expect_equal(as.vector(fit_large$coefficients$gamma), as.vector(fit_infinite$coefficients$gamma))
  expect_equal(as.vector(fit_large$coefficients$beta0), as.vector(fit_infinite$coefficients$beta0))
})

test_that("gibbsm works with infinite saturation for Exponential/Geyer potentials", {
  ntypes <- 4
  short_range <- matrix(0.01, ntypes, ntypes)
  medium_range <- matrix(0.02, ntypes, ntypes)
  long_range <- matrix(0.04, ntypes, ntypes)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 1e6

  set.seed(1)
  z <- ppjsdm::rppp(lambda = rep(1e2, ntypes))

  set.seed(1)
  fit_large <- gibbsm(z, short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      use_glmnet = FALSE,
                      max_dummy = max_dummy,
                      dummy_factor = dummy_factor,
                      saturation = 1e6,
                      model = "exponential",
                      medium_range_model = "Geyer")

  set.seed(1)
  fit_infinite <- gibbsm(z, short_range = short_range,
                         medium_range = medium_range,
                         long_range = long_range,
                         use_glmnet = FALSE,
                         max_dummy = max_dummy,
                         dummy_factor = dummy_factor,
                         saturation = Inf,
                         model = "exponential",
                         medium_range_model = "Geyer")

  expect_equal(as.vector(fit_large$coefficients$alpha), as.vector(fit_infinite$coefficients$alpha))
  expect_equal(as.vector(fit_large$coefficients$gamma), as.vector(fit_infinite$coefficients$gamma))
  expect_equal(as.vector(fit_large$coefficients$beta0), as.vector(fit_infinite$coefficients$beta0))
})

test_that("gibbsm works with infinite saturation for Exponential/Half exponential potentials", {
  ntypes <- 4
  short_range <- matrix(0.01, ntypes, ntypes)
  medium_range <- matrix(0.02, ntypes, ntypes)
  long_range <- matrix(0.04, ntypes, ntypes)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 1e6

  set.seed(1)
  z <- ppjsdm::rppp(lambda = rep(1e2, ntypes))

  set.seed(1)
  fit_large <- gibbsm(z, short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      use_glmnet = FALSE,
                      max_dummy = max_dummy,
                      dummy_factor = dummy_factor,
                      saturation = 1e6,
                      model = "exponential",
                      medium_range_model = "half_exponential")

  set.seed(1)
  fit_infinite <- gibbsm(z, short_range = short_range,
                         medium_range = medium_range,
                         long_range = long_range,
                         use_glmnet = FALSE,
                         max_dummy = max_dummy,
                         dummy_factor = dummy_factor,
                         saturation = Inf,
                         model = "exponential",
                         medium_range_model = "half_exponential")

  expect_equal(as.vector(fit_large$coefficients$alpha), as.vector(fit_infinite$coefficients$alpha))
  expect_equal(as.vector(fit_large$coefficients$gamma), as.vector(fit_infinite$coefficients$gamma))
  expect_equal(as.vector(fit_large$coefficients$beta0), as.vector(fit_infinite$coefficients$beta0))
})

test_that("gibbsm works with infinite saturation for Exponential/Square exponential potentials", {
  ntypes <- 4
  short_range <- matrix(0.01, ntypes, ntypes)
  medium_range <- matrix(0.02, ntypes, ntypes)
  long_range <- matrix(0.04, ntypes, ntypes)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 1e6

  set.seed(1)
  z <- ppjsdm::rppp(lambda = rep(1e2, ntypes))

  set.seed(1)
  fit_large <- gibbsm(z, short_range = short_range,
                      medium_range = medium_range,
                      long_range = long_range,
                      use_glmnet = FALSE,
                      max_dummy = max_dummy,
                      dummy_factor = dummy_factor,
                      saturation = 1e6,
                      model = "exponential",
                      medium_range_model = "square_exponential")

  set.seed(1)
  fit_infinite <- gibbsm(z, short_range = short_range,
                         medium_range = medium_range,
                         long_range = long_range,
                         use_glmnet = FALSE,
                         max_dummy = max_dummy,
                         dummy_factor = dummy_factor,
                         saturation = Inf,
                         model = "exponential",
                         medium_range_model = "square_exponential")

  expect_equal(as.vector(fit_large$coefficients$alpha), as.vector(fit_infinite$coefficients$alpha))
  expect_equal(as.vector(fit_large$coefficients$gamma), as.vector(fit_infinite$coefficients$gamma))
  expect_equal(as.vector(fit_large$coefficients$beta0), as.vector(fit_infinite$coefficients$beta0))
})
