#' Return a configuration representing the locations of 'ten' species of trees in a given CSIRO experimental plot index by 'index', in a given 'year'.
#'
#' @param index An index between 1 and 20 representing the experimental plot.
#' @param year A number between 1971 and 2013 representing the year.
#' @param prevalent Number of most prevalent species to use.
#' @param jitter Variance of the shifts applied to the locations of the trees.
#' @param average_height Return trait corresponding to average heights?
#' @param average_dbh Return trait corresponding to average stem diameters?
#' @importFrom stats rnorm
#' @export
get_qld_trees_configuration <- function(index = 15, year = 2013, prevalent = 10, jitter = 0.3, average_height = FALSE, average_dbh = FALSE) {
  stopifnot(index >= 1, index <= 20, year >= 1971, year <= 2013)
  raw_qld_trees <- ppjsdm::raw_qld_trees

  index_epNumber <- unique(raw_qld_trees$epNumber)[index]
  message("The chosen index corresponds to ", paste(index_epNumber, collapse = ", "), ".")
  modified_qld_trees <- raw_qld_trees[raw_qld_trees$epNumber %in% index_epNumber, ]

  used_year <- year
  repeat {
    if(year < 1971) {
      stop("There was no data for any year less than or equal to the requested year.")
    }
    condition <- modified_qld_trees$year == used_year
    if(any(condition)) {
      modified_qld_trees <- modified_qld_trees[condition, ]
      break
    } else {
      used_year <- used_year - 1
    }
  }

  if(used_year != year) {
    warning("There was no data for the requested year ", year, ", using year ", used_year, " instead.")
  }


  sp_prevalence <- rev(sort(table(modified_qld_trees$taxon)))
  sp_keep <- names(sp_prevalence)[1:prevalent]
  modified_qld_trees <- modified_qld_trees[modified_qld_trees$taxon %in% sp_keep, ]

  number_locations <- length(modified_qld_trees$coordinates_x_metres)
  jiggle <- rnorm(2 * number_locations, 0, jitter)
  modified_qld_trees$coordinates_x_metres <- modified_qld_trees$coordinates_x_metres + jiggle[1:number_locations]
  modified_qld_trees$coordinates_y_metres <- modified_qld_trees$coordinates_y_metres + jiggle[-(1:number_locations)]

  # Truncate everything to the observation region (it could have moved outside because of the jitter).
  modified_qld_trees$coordinates_x_metres[modified_qld_trees$coordinates_x_metres > 100] <- 100
  modified_qld_trees$coordinates_x_metres[modified_qld_trees$coordinates_x_metres < 0] <- 0
  modified_qld_trees$coordinates_y_metres[modified_qld_trees$coordinates_y_metres > 50] <- 50
  modified_qld_trees$coordinates_y_metres[modified_qld_trees$coordinates_y_metres < 0] <- 0

  taxon_factor <- factor(modified_qld_trees$taxon)
  configuration <- Configuration(modified_qld_trees$coordinates_x_metres, modified_qld_trees$coordinates_y_metres, taxon_factor)
  if(average_height) {
    v <- sapply(1:prevalent, function(i) mean(modified_qld_trees$height_metres[taxon_factor == levels(taxon_factor)[i]], na.rm = TRUE), USE.NAMES = FALSE)
    average_height_trait <- diag(v)
  }
  if(average_dbh) {
    v <- sapply(1:prevalent, function(i) mean(modified_qld_trees$dbh_centimetres[taxon_factor == levels(taxon_factor)[i]], na.rm = TRUE), USE.NAMES = FALSE)
    average_dbh_trait <- diag(v)
  }
  # TODO: I'm sure there's a better way to do this...
  if(average_height) {
    if(average_dbh) {
      list(configuration = configuration, average_height = average_height_trait, average_dbh = average_dbh_trait)
    } else {
      list(configuration = configuration, average_height = average_height_trait)
    }
  } else {
    if(average_dbh) {
      list(configuration = configuration, average_dbh = average_dbh_trait)
    } else {
      configuration
    }
  }
}
