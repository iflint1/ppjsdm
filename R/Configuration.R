#' Configuration constructor.
#'
#' @param x Values along the x-axis.
#' @param y Values along the y-axis.
#' @param types Types of the points.
#' @param marks Marks of the points.
#' @importFrom methods is
#' @export
Configuration <- local({
  function(x, y, types, marks) {
    if(nargs() == 1 && !missing(x)) {
      configuration <- as.Configuration(x)
    } else {
      if(!missing(x)) {
        n <- length(x)
      } else if(!missing(y)) {
        n <- length(y)
      } else if(!missing(types)) {
        n <- length(types)
      } else if(!missing(marks)) {
        n <- length(marks)
      }
      if(missing(x)) {
        x <- rep(1, n)
      }
      if(missing(y)) {
        y <- rep(1, n)
      }
      if(length(x) != length(y) || !is.vector(x) || !is.vector(y) || !is.numeric(x) || !is.numeric(y)) {
        stop("x and y should be numeric vectors with the same length.")
      }
      if(missing(types)) {
        types <- factor(rep("default", n))
      } else if(!is.factor(types)) {
        if(length(types) == n) {
          types <- factor(types)
        } else {
          types <- factor(rep(types, n))
        }
      } else if(length(types) != n) {
        stop("The types argument should be a factor with the same length as that of the other parameters.")
      }

      if(missing(marks)) {
        marks <- rep(1.0, n)
      } else if(length(marks) != n || !is.vector(marks) || !is.numeric(marks)) {
        print(marks)
        stop("The marks argument should be a numeric vector with the same length as that of the other parameters.")
      }
      configuration <- structure(list(x = x, y = y, types = types, marks = marks), class = "Configuration")
    }
    if(has_duplicates(configuration)) {
      warning("There are duplicate points in the configuration.")
    }
    configuration
  }
})

#' Convert to our configuration class.
#'
#' @param configuration Configuration.
#' @export
as.Configuration <- function(configuration) {
  UseMethod("as.Configuration", configuration)
}

#' @method as.Configuration Configuration
as.Configuration.Configuration <- function(configuration) {
  configuration
}

#' @method as.Configuration ppp
#' @export
as.Configuration.ppp <- function(configuration) {
  x <- configuration$x
  y <- configuration$y
  marks <- configuration$marks
  if(is.null(marks)) {
    types <- factor(rep("default", length(configuration$x)))
  } else {
    types <- marks
  }
  Configuration(x = x,
                y = y,
                types = types,
                marks = rep(1.0, length(x)))
}

#' @method as.Configuration numeric
as.Configuration.numeric <- function(configuration) {
  if(is.vector(configuration)) {
    x <- configuration
    Configuration(x = x,
                  y = rep(1, length(x)),
                  types = factor(rep("default", length(x))),
                  marks = rep(1, length(x)))
  } else {
    as.Configuration.default(configuration)
  }
}

#' @method as.Configuration matrix
as.Configuration.matrix <- function(configuration) {
  if(ncol(configuration) < 2) {
    stop("Matrix cannot be converted to a configuration.")
  }
  x <- configuration[, 1]
  y <- configuration[, 2]
  if(ncol(configuration) == 2) {
    marks <- rep(1, length(x))
  } else {
    marks <- configuration[, 3]
  }
  Configuration(x = x,
                y = y,
                types = factor(rep("default", length(x))),
                marks = marks)
}

#' @method as.Configuration default
as.Configuration.default <- function(configuration) {
  configuration <- as.data.frame(configuration)
  x_indices <- which(names(configuration) == "x")
  y_indices <- which(names(configuration) == "y")
  if(length(x_indices) == 0) {
    if(length(y_indices) == 0) {
      stop("Cannot convert to configuration.")
    } else {
      y <- configuration[, y_indices]
      x <- rep(1, length(y))
    }
  } else {
    x <- configuration[, x_indices]
    if(length(y_indices) == 0) {
      y <- rep(1, length(x))
    } else {
      y <- configuration[, y_indices]
    }
  }
  types_indices <- which(names(configuration) == "types")
  marks_indices <- which(names(configuration) == "marks")
  if(length(types_indices) == 0) {
    types <- factor(rep("default", length(x)))
  } else {
    types <- configuration[, types_indices]
  }
  if(length(marks_indices) == 0) {
    marks <- rep(1, length(x))
  } else {
    marks <- configuration[, marks_indices]
  }
  Configuration(x = x,
                y = y,
                types = types,
                marks = marks)
}

#' Subset of the configuration consisting only of points of the i-th type.
#'
#' @param x Configuration.
#' @param i Subscript index.
#' @export
"[.Configuration" <- function(x, i) {
  types <- x$types
  subset_indices <- types %in% levels(types)[as.integer(i)]
  Configuration(x = x$x[subset_indices], y = x$y[subset_indices], types = factor(types[subset_indices]), marks = x$marks[subset_indices])
}

format <- function(configuration) {
  number_points <- length(configuration$x)
  str <- paste0("An S3 object representing a configuration.\n\nNumber of points: ",
                paste0(number_points, collapse = ", "),
                ".\n")
  if(number_points > 0 && number_points < 50) {
    str <- paste0(str, "\nPoints in the format (x-coordinate, y-coordinate, type): ",
                  paste0("(", configuration$x, ", ", configuration$y, ", ", configuration$types, ")", collapse = ", "),
                  ".\n\nCoordinate(s) along the x-axis: ",
                  paste0(configuration$x, collapse = ", "),
                  ".\n\nCoordinate(s) along the y-axis: ",
                  paste0(configuration$y, collapse = ", "),
                  ".\n\nType(s) of the point(s): ",
                  paste0(configuration$types, collapse = ", "),
                  ".\n\nMarks(s) of the point(s): ",
                  paste0(configuration$marks, collapse = ", "),
                  ".\n")
  }
  str
}

#' @method print Configuration
#' @export
print.Configuration <- function(x, ...) {
  cat(format(x))
}

#' @importFrom graphics legend par plot
#' @method plot Configuration
#' @export
plot.Configuration <- function(x, window, ...) {
  if(length(x$x) > 0) {
    if(missing(window)) {
      x_range <- c(min(x$x), max(x$x))
      y_range <- c(min(x$y), max(x$y))
    } else {
      x_range <- x_range(window)
      y_range <- y_range(window)
    }
    # TODO: Allow for user to supply colors & labels.
    plot(x$x,
         x$y,
         xlim = x_range,
         ylim = y_range,
         col = droplevels(x$types),
         xlab = "x",
         ylab = "y",
         main = "Points in the configuration",
         ... )

    par(mar = c(5, 4, 4, 20) + 0.1)

    legend("topleft",
           xpd = TRUE,
           inset = c(1.03, 0),
           legend = levels(droplevels(x$types)),
           col = 1:length(levels(x$types)),
           pch = 1,
           pt.cex = 1,
           cex = 0.8,
           title = "Types of the points")
  }
}

#' Access x-coordinates of a configuration
#'
#' @param configuration The configuration.
#' @export
x_coordinates <- function(configuration) {
  configuration$x
}

#' Access y-coordinates of a configuration
#'
#' @param configuration The configuration.
#' @export
y_coordinates <- function(configuration) {
  configuration$y
}

#' Access types of a configuration
#'
#' @param configuration The configuration.
#' @export
types <- function(configuration) {
  configuration$types
}

#' Access marks of a configuration
#'
#' @param configuration The configuration.
#' @export
marks <- function(configuration) {
  configuration$marks
}

#' Convert a configuration class to a ppp from the SpatStat package.
#'
#' @param X Configuration.
#' @param W Window on which the points are located.
#' @param ... Not used.
#' @param fatal Not used.
#' @importFrom spatstat as.owin owin ppp
#' @export
as.ppp.Configuration <- function(X, W, ..., fatal = TRUE) {
  types <- X$types
  # The lines below "unfactor" types
  types <- as.character(types)
  if(all(sapply(types, function(t) suppressWarnings(!is.na(as.numeric(t)))))) {
    types <- as.numeric(types)
  }
  ppp(X$x, X$y, window = as.owin(W), marks = types)
}

#' Number of points in a configuration
#'
#' @param x Configuration.
#' @export
length.Configuration <- function(x) {
  length(x$x)
}
