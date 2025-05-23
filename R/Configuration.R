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
    if(nargs() == 0) {
      configuration <- as.Configuration(vector(mode = "numeric", length = 0))
    } else if(nargs() == 1 && !missing(x)) {
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
  if(is.character(i)) {
    target_types <- i
  } else {
    target_types <- levels(types)[as.integer(i)]
  }
  subset_indices <- types %in% target_types
  Configuration(x = x$x[subset_indices],
                y = x$y[subset_indices],
                types = factor(types[subset_indices]),
                marks = x$marks[subset_indices])
}

format_configuration <- function(configuration) {
  number_points <- length(configuration$x)
  min_marks <- if(number_points > 0) {
    min(configuration$marks)
  } else {
    -Inf
  }
  max_marks <- if(number_points > 0) {
    max(configuration$marks)
  } else {
    Inf
  }
  str <- paste0("A configuration of points.\n\nTotal number of points: ",
                number_points,
                ".\nTypes: ",
                paste0(levels(configuration$types), collapse = ", "),
                ".\nMarks in range [",
                min_marks,
                ", ",
                max_marks,
                "].\nNumber of points by type: ",
                paste0(paste0(levels(configuration$types), " = ", table(configuration$types)), collapse = ", "),
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
  cat(format_configuration(x))
}

#' @importFrom ggplot2 aes coord_equal element_text geom_point ggplot ggtitle guides scale_colour_manual scale_colour_viridis_d scale_shape_manual theme_minimal xlab xlim ylab ylim
#' @importFrom rlang .data
#' @importFrom scales breaks_extended
#' @method plot Configuration
#' @export
plot.Configuration <- function(x,
                               window,
                               colours,
                               shapes,
                               base_size = 12,
                               types_name = "Types",
                               mark_range = c(1, 6), ...) {
  if(length(x$x) > 0) {
    if(missing(window)) {
      x_range <- c(min(x$x), max(x$x))
      y_range <- c(min(x$y), max(x$y))
    } else {
      window <- as.Window(window)
      x_range <- x_range(window)
      y_range <- y_range(window)
    }

    df <- data.frame(x = x$x, y = x$y, Marks = x$marks)

    df[[types_name]] <- x$types

    if(missing(shapes)) {
      shapes <- rep(c(16, 17, 15, 18), length.out = nlevels(df[[types_name]]))
    }

    g <- ggplot(data = df)

    if(!all(df$Marks == 1.)) {
      g <- g + geom_point(aes(x = .data$x, y = .data$y, colour = .data[[types_name]],
                              shape = .data[[types_name]], size = .data$Marks), alpha = 0.8) +
        scale_size(range = mark_range, breaks = breaks_extended(6))
      nr <- 6
    } else {
      g <- g + geom_point(aes(x = .data$x, y = .data$y, colour = .data[[types_name]],
                              shape = .data[[types_name]]), size = 2.5, alpha = 0.8)
      nr <- 8
    }

    g <- g + xlim(x_range[1], x_range[2]) +
      ylim(y_range[1], y_range[2]) +
      scale_shape_manual(values = shapes) +
      xlab(NULL) +
      ylab(NULL) +
      ggtitle("") +
      coord_equal() +
      theme_minimal(base_size = base_size)

    if(!missing(colours)) {
      g <- g + scale_colour_manual(values = rep(colours, nlevels(df[[types_name]])))
    } else {
      g <- g + scale_colour_viridis_d(end = 0.9, option = "turbo")
    }

    g <- g + guides(colour = guide_legend(order = 1,
                                          nrow = nr,
                                          override.aes = list(size = 5)),
                    shape = guide_legend(order = 1,
                                         nrow = nr),
                    size =  guide_legend(order = 2))

    g
  }
}

#' Access x-coordinates of a configuration
#'
#' @param configuration Configuration.
#' @export
x_coordinates <- function(configuration) {
  configuration$x
}

#' Access y-coordinates of a configuration
#'
#' @param configuration Configuration.
#' @export
y_coordinates <- function(configuration) {
  configuration$y
}

#' Access types of a configuration
#'
#' @param configuration Configuration.
#' @export
types <- function(configuration) {
  configuration$types
}

#' Access marks of a configuration
#'
#' @param x Configuration.
#' @param ... Unused.
#' @importFrom spatstat.geom marks
#' @exportS3Method spatstat.geom::marks Configuration
#' @examples
#' set.seed(1)
#'
#' # Create a configuration
#' configuration <- ppjsdm::Configuration(x = 1:4, y = 2:5, marks = runif(4))
#'
#' # Get its marks
#' print(marks(configuration))
#'
marks.Configuration <- function(x, ...) {
  x$marks
}

#' Convert Configuration to a data.frame.
#'
#' @param x Configuration.
#' @param ... Unused.
#' @method as.data.frame Configuration
#' @examples
#' set.seed(1)
#'
#' # Create a configuration
#' configuration <- ppjsdm::Configuration(x = 1:4, y = 2:5, marks = runif(4))
#'
#' # Convert it to a data.frame
#' print(as.data.frame(configuration))
#'
#' @export
as.data.frame.Configuration <- function(x, ...) {
  data.frame(x = x$x,
             y = x$y,
             types = x$types,
             marks = x$marks)
}

#' Convert a configuration class to a ppp from the SpatStat package.
#'
#' @param X Configuration.
#' @param W Window on which the points are located.
#' @param ... Not used.
#' @param fatal Not used.
#' @importFrom spatstat.geom as.owin owin ppp
#' @export
as.ppp.Configuration <- function(X, W, ..., fatal = TRUE) {
  ppp(x = X$x, y = X$y, window = as.owin(W), marks = X$types)
}

#' Number of points in a configuration
#'
#' @param x Configuration.
#' @method length Configuration
#' @export
length.Configuration <- function(x) {
  length(x$x)
}
