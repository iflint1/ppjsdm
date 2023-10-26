#' Area of a Window
#'
#' Imported directly from `spatstat`, in order to provide methods for our window classes,
#' see for example `?area.Rectangle_window`. The documentation here is sparse, and we
#' refer to the `spatstat` documentation for more details.
#' @param w Window.
#' @importFrom spatstat.geom area
#' @export
#' @md
area <- area

#' Volume of an Object
#'
#' Imported directly from `spatstat`, in order to provide methods for our window classes,
#' see for example `?volume.Rectangle_window`. The documentation here is sparse, and we
#' refer to the `spatstat` documentation for more details.
#' @param x Window.
#' @importFrom spatstat.geom volume
#' @export
#' @md
volume <- volume

#' Convert Data To Class owin
#'
#' Imported directly from `spatstat`, in order to provide methods for our window classes,
#' see for example `?as.owin.Rectangle_window`. The documentation here is sparse, and we
#' refer to the `spatstat` documentation for more details.
#' @param W Window.
#' @param ... Ignored.
#' @param fatal What to do if the data cannot be converted to an observation window?
#' @importFrom spatstat.geom as.owin
#' @export
#' @md
as.owin <- as.owin

#' Marks of a Point Pattern
#'
#' Imported directly from `spatstat`, in order to provide a method for our Configuration class, see
#' `?marks.Configuration`. The documentation here is sparse, and we
#' refer to the `spatstat` documentation for more details.
#' @param x Configuration of points.
#' @param ... Ignored.
#' @importFrom spatstat.geom marks
#' @export
#' @md
marks <- marks
