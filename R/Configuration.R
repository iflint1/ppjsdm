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
      if(length(x) != length(y) || !is.vector(x) || !is.vector(y) || !is.numeric(x) || !is.numeric(y)) {
        stop("x and y should be numeric vectors with the same length.")
      }
      if(missing(types)) {
        types <- factor(rep("default", length(x)))
      } else if(!is.factor(types)) {
        if(length(types) == length(x)) {
          types <- factor(types)
        } else {
          types <- factor(rep(types, length(x)))
        }
      }
      if(length(types) != length(x)) {
        stop("The types argument should be a factor with the same length as that of x.")
      }

      if(missing(marks)) {
        marks <- rep(1.0, length(x))
      } else {
        if(length(marks) != length(x) || !is.vector(marks) || !is.numeric(marks)) {
          stop("The marks argument should be a numeric vector with the same length as that of x.")
        }
      }
      configuration <- structure(list(x = x, y = y, types = types, marks = marks), class = "Configuration")
    }
    if(has_duplicates(configuration)) {
      warning("There are duplicate points in the configuration.")
    }
    configuration
  }
})


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
  str <- paste0("An S3 object representing a configuration.\n\n")
  str <- paste0(str, "Number of points: ")
  str <- paste0(str, paste0(number_points, collapse = ", "))
  str <- paste0(str, ".\n")
  if(number_points > 0 && number_points < 50) {
    str <- paste0(str, "\nPoints in the format (x-coordinate, y-coordinate, type): ")
    str <- paste0(str, paste0("(", configuration$x, ", ", configuration$y, ", ", configuration$types, ")", collapse = ", "))
    str <- paste0(str, ".\n\nCoordinate(s) along the x-axis: ")
    str <- paste0(str, paste0(configuration$x, collapse = ", "))
    str <- paste0(str, ".\n\nCoordinate(s) along the y-axis: ")
    str <- paste0(str, paste0(configuration$y, collapse = ", "))
    str <- paste0(str, ".\n\nType(s) of the point(s): ")
    str <- paste0(str, paste0(configuration$types, collapse = ", "))
    str <- paste0(str, ".\n\nMarks(s) of the point(s): ")
    str <- paste0(str, paste0(configuration$marks, collapse = ", "))
    str <- paste0(str, ".\n")
  }
  str
}

#' @method print Configuration
print.Configuration <- function(x, ...) {
  str <- format(x)
  cat(str)
}

#' @importFrom graphics legend par plot
#' @method plot Configuration
plot.Configuration <- function(x, window = NULL, ...) {
  if(length(x$x) > 0) {
    if(missing(window)) {
      x_range <- c(min(x_coordinates(x)), max(x_coordinates(x)))
      y_range <- c(min(y_coordinates(x)), max(y_coordinates(x)))
    } else {
      x_range <- x_range(window)
      y_range <- y_range(window)
    }
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
  ppp(X$x, X$y, window = as.owin(W), marks = X$types)
}

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
as.Configuration.ppp <- function(configuration) {
  marks <- configuration$marks
  if(is.null(marks)) {
    Configuration(x = configuration$x,
                  y = configuration$y,
                  types = factor(rep("default", length(configuration$x))),
                  marks = rep(1, length(configuration$x)))
  } else {
    Configuration(x = configuration$x,
                  y = configuration$y,
                  types = marks,
                  marks = rep(1, length(configuration$x)))
  }
}

#' @method as.Configuration vector
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

#' @method as.Configuration default
as.Configuration.default <- function(configuration) {
  Configuration(x = configuration[, 1],
                y = configuration[, 2],
                types = factor(rep("default", length(configuration[, 1]))),
                marks = rep(1, length(configuration[, 1])))
}


#' Add a point to a configuration.
#'
#' @param configuration Configuration.
#' @param x Point.
#' @param type Point type.
#' @param mark Point mark
#' @export
add_to_configuration <- function(configuration, x, type = "default", mark = 1) {
  if(!is.vector(x) | !is.numeric(x) | length(x) != 2) {
    stop("Wrong format for x.")
  }
  configuration$x <- c(configuration$x, x[1])
  configuration$y <- c(configuration$y, x[2])
  configuration$types <- unlist(list(configuration$types, factor(type)))
  configuration$marks <- mark
  configuration
}

#' Remove a point from a configuration by index.
#'
#' @param configuration Configuration.
#' @param index Index.
#' @export
remove_from_configuration_by_index <- function(configuration, index) {
  configuration$x <- configuration$x[-index]
  configuration$y <- configuration$y[-index]
  configuration$types <- droplevels(configuration$types[-index])
  configuration$marks <- configuration$marks[-index]
  configuration
}


#' Remove a point from a configuration (if it is in the configuration).
#'
#' @param configuration Configuration.
#' @param x Point.
#' @param type Point type.
#' @param mark Point mark.
#' @export
remove_from_configuration <- function(configuration, x, type = "default", mark = 1) {
  index <- match(x[1], configuration$x)
  if(!is.na(index) & configuration$y[index] == x[2] & configuration$types[index] == type & configuration$marks[index] == mark) {
    remove_from_configuration_by_index(configuration, index)
  } else {
    configuration
  }
}

#' Number of points in a configuration
#'
#' @param x Configuration.
#' @export
length.Configuration <- function(x) {
  length(x$x)
}

#' Check whether a configuration is empty.
#'
#' @param configuration Configuration.
#' @export
is_empty <- function(configuration) {
  length(configuration) == 0
}

#' Remove random point from the configuration.
#'
#' @param configuration Configuration.
#' @export
remove_random_point <- function(configuration) {
  index <- sample(length(configuration$x), 1)
  location <- c(configuration$x[index], configuration@y[index])
  type <- configuration$types[index]
  mark <- configuration$marks[index]
  configuration$x <- configuration$x[-index]
  configuration$y <- configuration$y[-index]
  configuration$types <- droplevels(configuration$types[-index])
  configuration$marks <- droplevels(configuration$marks[-index])
  list(configuration = configuration, location = location, type = type, mark = mark)
}
