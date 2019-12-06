#' Return a configuration representing the locations of 'ten' species of trees in a given CSIRO experimental plot index by 'index', in a given 'year'.
#'
#' @param index An index between 1 and 20 representing the experimental plot.
#' @param year A number between 1971 and 2013 representing the year.
#' @param prevalent Number of most prevalent species to use.
#' @param jitter Variance of the shifts applied to the locations of the trees.
#' @export
get_csiro_trees_configuration <- function(index = 15, year = 2013, prevalent = 10, jitter = 0.3) {
  stopifnot(index >= 1, index <= 20, year >= 1971, year <= 2013)
  raw_csiro_trees <- ppjsdm::raw_csiro_trees

  index_epNumber <- unique(raw_csiro_trees$epNumber)[index]
  message("The chosen index corresponds to ", paste(index_epNumber, collapse = ", "), ".")
  modified_csiro_trees <- raw_csiro_trees[raw_csiro_trees$epNumber %in% index_epNumber, ]

  used_year <- year
  repeat {
    if(year < 1971) {
      stop("There was no data for any year less than or equal to the requested year.")
    }
    condition <- modified_csiro_trees$year == used_year
    if(any(condition)) {
      modified_csiro_trees <- modified_csiro_trees[condition, ]
      break
    } else {
      used_year <- used_year - 1
    }
  }

  if(used_year != year) {
    warning("There was no data for the requested year ", year, ", using year ", used_year, " instead.")
  }


  sp_prevalence <- rev(sort(table(modified_csiro_trees$taxon)))
  sp_keep <- names(sp_prevalence)[1:prevalent]
  modified_csiro_trees <- modified_csiro_trees[modified_csiro_trees$taxon %in% sp_keep, ]

  number_locations <- length(modified_csiro_trees$coordinates_x_metres)
  jiggle <- rnorm(2 * number_locations, 0, jitter)
  modified_csiro_trees$coordinates_x_metres <- modified_csiro_trees$coordinates_x_metres + jiggle[1:number_locations]
  modified_csiro_trees$coordinates_y_metres <- modified_csiro_trees$coordinates_y_metres + jiggle[-(1:number_locations)]

  Configuration(modified_csiro_trees$coordinates_x_metres, modified_csiro_trees$coordinates_y_metres, factor(modified_csiro_trees$taxon))
}
