# TODO: Test factors

#' A class representing a configuration.
#'
#' @param x Values along the x-axis.
#' @param y Values along the y-axis.
#' @param types Types of the points.
#' @export
setClass(Class = "Configuration", slots = list(x = "numeric", y = "numeric", types = "factor"))


#' Configuration constructor.
#'
#' @param x Values along the x-axis.
#' @param y Values along the y-axis.
#' @param types Types of the points.
#' @importFrom methods is new
#' @export
Configuration <- local({
  function(x, y, types) {
    if(nargs() == 1 && !missing(x)) {
      if(is(x, "Configuration")) {
        configuration <- x
      } else {
        configuration <- new("Configuration", x = x[, 1], y = x[, 2], types = factor(rep("default", length(x[, 1]))))
      }
    } else {
      if(length(x) != length(y) || !is.vector(x) || !is.vector(y)) {
        stop("x and y should be vectors with the same length.")
      }
      if(missing(types)) {
        types <- factor(rep("default", length(x)))
      }
      configuration <- new("Configuration", x = x, y = y, types = types)
    }
    if(has_duplicates(configuration)) {
      warning("There are duplicate points in the configuration.")
    }
    configuration
  }
})


#' Configuration subset operator.
#'
#' @param x Configuration.
#' @param i Subscript index.
#' @param j Ignored.
#' @param ... Ignored.
#' @param drop Ignored.
#' @importFrom methods initialize
#' @export
setMethod("[", c("Configuration", "numeric", "missing", "ANY"),
  function(x, i, j, ..., drop = TRUE) {
    types <- x@types
    subset_indices <- types %in% levels(types)[as.integer(i)]
    initialize(x, x = x@x[subset_indices], y = x@y[subset_indices], types = factor(types[subset_indices]))
  }
)

#' Print a configuration class.
#' @param object Configuration to print.
#' @importFrom methods show
#' @export
setMethod("show", "Configuration", function(object) {
  number_points <- length(object@x)
  cat("An S4 object representing a configuration.\n\n")
  cat("Number of points:", paste0(number_points, ".\n"))
  if(number_points > 0 && number_points < 50) {
    cat("\n")
    cat("Points in the format (x-coordinate, y-coordinate, type): ")
    cat(paste0("(", object@x, ", ", object@y, ", ", object@types, ")"), sep =", ")
    cat(".\n\n")
    cat("Coordinate(s) along the x-axis: ")
    cat(object@x, sep = ", ")
    cat(".\n\n")
    cat("Coordinate(s) along the y-axis: ")
    cat(object@y, sep = ", ")
    cat(".\n\n")
    cat("Type(s) of the point(s): ")
    cat(object@types, sep = ", ")
    cat(".\n")
  }
})

#' Number of points in a configuration
#' @param configuration The configuration.
#' @export
get_number_points <- function(configuration) {
  # TODO: The code below is horrible, it can definitely be vectorised...
  types <- configuration@types
  number_types <- length(levels(types))
  result <- rep(NA, number_types)
  index <- 1
  for(i in levels(types)) {
    result[index] <- length(configuration@x[configuration@types == i])
    index <- index + 1
  }
  result
}


#' Plot a configuration class.
#' @param x Configuration to plot.
#' @param window Window the configuration belongs to.
#' @param ... Other arguments not yet used.
#' @importFrom graphics plot par legend
#' @export
plot.Configuration <- function(x, window = NULL, ...) {
  if(sum(get_number_points(x)) > 0) {
    if(missing(window)) {
      plot(x@x, x@y, col = x@types, xlab = "x", ylab = "y", main = "Points in the configuration")
    } else {
      plot(x@x, x@y, xlim = window@x_range, ylim = window@y_range, col = x@types, xlab = "x", ylab = "y", main = "Points in the configuration")
    }
    # Add extra space to the right of plot area; change clipping to figure
    par(mar = c(5, 4, 4, 30) + 0.1, xpd = TRUE)

    legend("topleft", inset = c(1.03, 0), legend = unique(x@types), col = 1:length(x@types), pch = 1, title = "Types of the points")
  }
}

#' Access x-coordinates of a configuration
#' @param configuration The configuration.
#' @export
x_coordinates <- function(configuration) {
  configuration@x
}

#' Access y-coordinates of a configuration
#' @param configuration The configuration.
#' @export
y_coordinates <- function(configuration) {
  configuration@y
}

#' Access types of a configuration
#' @param configuration The configuration.
#' @export
types <- function(configuration) {
  configuration@types
}

#' Convert a configuration class to a ppp from the SpatStat package.
#' @param X Configuration.
#' @param W Window on which the points are located.
#' @param ... Not used.
#' @param fatal Not used.
#' @importFrom spatstat as.owin owin ppp
#' @export
as.ppp.Configuration <- function(X, W, ..., fatal = TRUE) {
  ppp(X@x, X@y, window = as.owin(W), marks = X@types)
}

#' Add a point to a configuration.
#' @param configuration Configuration.
#' @param x Point.
#' @param type Point type.
#' @export
add_to_configuration <- function(configuration, x, type) {
  configuration@x <- c(configuration@x, x[1])
  configuration@y <- c(configuration@y, x[2])
  configuration@types <- unlist(list(configuration@types, factor(type)))
  configuration
}

#' Remove a point from configuration (if it is in the configuration).
#' @param configuration Configuration.
#' @param x Point.
#' @param type Point type.
#' @export
remove_from_configuration <- function(configuration, x, type) {
  index <- match(x[1], configuration@x)
  if(!is.na(index) && configuration@y[index] == x[2] && configuration@types[index] == type) {
    configuration@x <- configuration@x[-index]
    configuration@y <- configuration@y[-index]
    configuration@types <- droplevels(configuration@types[-index])
  }
  configuration
}

#' Check whether a configuration is empty.
#' @param configuration Configuration.
#' @export
is_empty <- function(configuration) {
  length(configuration@x) == 0
}

#' Remove random point from the configuration.
#' @param configuration Configuration.
#' @export
remove_random_point <- function(configuration) {
  index <- sample(length(configuration@x), 1)
  location <- c(configuration@x[index], configuration@y[index])
  type <- configuration@types[index]
  configuration@x <- configuration@x[-index]
  configuration@y <- configuration@y[-index]
  configuration@types <- droplevels(configuration@types[-index])
  list(configuration = configuration, location = location, type = type)
}
