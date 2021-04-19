library(ppjsdm)

context("gibbsm")

test_that("gibbsm returns correct parameters for Geyer potentials potentials", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "Geyer",
                medium_range_model = "Geyer")

  target_alpha <- cbind(c(-6.59604203690729029574413289083167910575866699218750, -0.28730698002687304892788233701139688491821289062500),
                        c(-0.28730698002687304892788233701139688491821289062500, 0.30729800429036019382422750823025126010179519653320))
  target_gamma <- cbind(c(-0.05161243384419380436645141685403359588235616683960, 0.05582244876153419355091500619892030954360961914062),
                        c(0.05582244876153419355091500619892030954360961914062, -0.01029334248016372521661310202034655958414077758789))
  target_beta0 <- c(4.56015654138443782272815951728262007236480712890625, 4.68522188044612519775000691879540681838989257812500)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
})

test_that("gibbsm returns correct parameters for Exponential/Geyer potentials", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "exponential",
                medium_range_model = "Geyer")

  target_alpha <- cbind(c(0.37477043101832158145114703984290827065706253051758, -0.19498177190133467173716041997977299615740776062012),
                        c(-0.19498177190133467173716041997977299615740776062012, 0.27385642132631649037222132392344065010547637939453))
  target_gamma <- cbind(c(-0.08229360683943837129206144709314685314893722534180, 0.08428574178534278815710933940863469615578651428223),
                        c(0.08428574178534278815710933940863469615578651428223, -0.04164321456553284062085396044494700618088245391846))
  target_beta0 <- c(4.48193857807983242480531771434471011161804199218750, 4.67031013303837738703805371187627315521240234375000)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
})

test_that("gibbsm returns correct parameters for Exponential/Half exponential potentials", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

  # Exponential/Half exponential potentials
  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "exponential",
                medium_range_model = "half_exponential")

  target_alpha <- cbind(c(0.47672377817360295226833954984613228589296340942383, -0.13114953861490216691088050993130309507250785827637),
                        c(-0.13114953861490216691088050993130309507250785827637, 0.26278957802394509180032855510944500565528869628906))
  target_gamma <- cbind(c(-0.34849770430217219541191298048943281173706054687500, 0.04682457184185797455233668529217538889497518539429),
                        c(0.04682457184185797455233668529217538889497518539429, -0.03036531945194816464739240302606049226596951484680))
  target_beta0 <- c(4.71716761058562195074728151666931807994842529296875, 4.67488892714819215257193718571215867996215820312500)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
})

test_that("gibbsm returns correct parameters for Exponential/Square exponential potentials", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "exponential",
                medium_range_model = "square_exponential")

  target_alpha <- cbind(c(0.51714305198667565033332493840134702622890472412109, -0.28437429917999018647023490302672144025564193725586),
                        c(-0.28437429917999018647023490302672144025564193725586, 0.33298307882591560158402899105567485094070434570312))
  target_gamma <- cbind(c(-0.20478739895648501168068378319730982184410095214844, 0.14560338667980046722938425318716326728463172912598),
                        c(0.14560338667980046722938425318716326728463172912598, -0.09353107827259252160523317343177041038870811462402))
  target_beta0 <- c(4.49886964648366038943549938267096877098083496093750, 4.67027070309024860961244485224597156047821044921875)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
})

test_that("gibbsm works with infinite saturation for Geyer potentials potentials", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 1e6

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

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
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

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
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

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
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

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
