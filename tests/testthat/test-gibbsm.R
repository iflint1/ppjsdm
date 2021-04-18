library(ppjsdm)

context("gibbsm")

test_that("gibbsm returns correct parameters for given test values", {
  short_range <- matrix(0.01, 2, 2)
  medium_range <- matrix(0.02, 2, 2)
  long_range <- matrix(0.04, 2, 2)
  max_dummy <- 1e3
  dummy_factor <- 1e6
  saturation <- 2

  set.seed(1)
  z <- ppjsdm::rppp(lambda = c(1e2, 1e2))

  # Geyer potentials
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

  # Exponential/Geyer potentials
  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "exponential",
                medium_range_model = "Geyer")

  target_alpha <- cbind(c(0.35010323872594129035462628962704911828041076660156, -0.17759870456817056227905027299129869788885116577148),
                        c(-0.17759870456817056227905027299129869788885116577148, 0.35256504184130565970889392701792530715465545654297))
  target_gamma <- cbind(c(-0.10084714698976385283124557190603809431195259094238, 0.09702105429207094622334750511072343215346336364746),
                        c(0.09702105429207094622334750511072343215346336364746, -0.03847266422492109927411263470276026055216789245605))
  target_beta0 <- c(4.48595848094682470019733955268748104572296142578125, 4.63547533681566203966895045596174895763397216796875)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))

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

  target_alpha <- cbind(c(0.41961698469347419315766956060542725026607513427734, -0.07094778573517783459845276183841633610427379608154),
                        c(-0.07094778573517783459845276183841633610427379608154, 0.20364521431050097710624413593905046582221984863281))
  target_gamma <- cbind(c(-0.42604101967695279240544437016069423407316207885742, 0.05816250390604346676148850292520364746451377868652),
                        c(0.05816250390604346676148850292520364746451377868652, -0.03335581619183852764010111968673299998044967651367))
  target_beta0 <- c(4.77603317684604888881949591450393199920654296875000, 4.67087840259726849723165287286974489688873291015625)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))

  # Exponential/Square exponential potentials
  fit <- gibbsm(z, short_range = short_range,
                medium_range = medium_range,
                long_range = long_range,
                use_glmnet = FALSE,
                max_dummy = max_dummy,
                dummy_factor = dummy_factor,
                saturation = saturation,
                model = "exponential",
                medium_range_model = "square_exponential")

  target_alpha <- cbind(c(0.38066332862224716571120097796665504574775695800781, -0.26417831434202587725934563422924838960170745849609),
                        c(-0.26417831434202587725934563422924838960170745849609, 0.35588460295789992038351101655280217528343200683594))
  target_gamma <- cbind(c(-0.25978657093834128799514360252942424267530441284180, 0.11162640512788384039577493922479334287345409393311),
                        c(0.11162640512788384039577493922479334287345409393311, -0.07293830700391636112644988543252111412584781646729))
  target_beta0 <- c(4.58641744896081160476342120091430842876434326171875, 4.66904973103962372960040738689713180065155029296875)

  expect_equal(as.vector(fit$coefficients$alpha), as.vector(target_alpha))
  expect_equal(as.vector(fit$coefficients$gamma), as.vector(target_gamma))
  expect_equal(as.vector(fit$coefficients$beta0), as.vector(target_beta0))
})
