#' Summary of a Fitted Gibbs point process.
#'
#' @param object A fitted model object.
#' @param ... Ignored.
#' @importFrom stats pnorm qnorm
#' @export
summary.gibbsm <- function(object, ...) {
  y <- list()
  class(y) <- "summary.gibbsm"

  vc <- vcov.gibbsm(object, ...)
  se <- sqrt(diag(vc))

  coefficients <- object$coefficients
  # TODO: Don't want to be doing the line below
  coefficients <- c(beta0 = coefficients$beta0, alpha = coefficients$alpha)
  one_ninetysix <- qnorm(0.975)
  lo <- coefficients - one_ninetysix * se
  hi <- coefficients + one_ninetysix * se
  zval <- coefficients / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  psig <- cut(pval,
              breaks = c(0,0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "  "),
              include.lowest = TRUE)
  y$coefficients <- data.frame(coefficients = coefficients,
                               se = se,
                               CI95_lo = lo,
                               CI95_hi = hi,
                               Ztest = psig,
                               Pval = pval,
                               Zval = zval)
  y
}

#' Print summary of a Fitted Gibbs point process.
#'
#' @param x A summary of a fitted model object.
#' @param ... Ignored.
#' @export
print.summary.gibbsm <- function(x, ...) {
  print(x$coefficients)
}
