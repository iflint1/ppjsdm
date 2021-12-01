#' Execute a function on a list/vector of objects, using educated guesses to execute as many calls
#' as possible within a specified time limit. The function is always called at least once (on the
#' first object).
#'
#' @param objects List/vector of objects.
#' @param func Function to apply to each object.
#' @param time_limit Time limit in `unit` that can be spent running this function.
#' @param unit Unit used to measure the time limit (hours, mins, secs, etc).
#'
#' @importFrom stats sd
#' @keywords internal
execute_until_time_limit <- function(objects, func, time_limit = Inf, unit = "hours") {
  objects_type <- class(objects)
  if(!is.list(objects) & !is.vector(objects)) {
    stop("objects' class not recognised.")
  }
  if(is.infinite(time_limit)) {
    if(is.list(objects)) {
      lapply(objects, func)
    } else {
      sapply(objects, func)
    }
  } else {
    result <- vector(mode = objects_type, length = length(objects))
    execution_times <- c()
    time_left <- time_limit
    for(i in seq_len(length(objects))) {
      if(length(execution_times) > 0) {
        if(length(execution_times) > 1) {
          max_time <- mean(execution_times) + 2 * sd(execution_times)
        } else { # execution_times only has one element
          max_time <- 2 * execution_times[1]
        }
      } else {
        max_time <- -Inf
      }
      if(max_time < time_left) {
        start <- Sys.time()
        if(is.list(objects)) {
          result[[i]] <- func(objects[[i]])
        } else {
          result[i] <- func(objects[i])
        }
        execution_time <- as.numeric(Sys.time() - start, unit = unit)
        execution_times <- c(execution_times, execution_time)
        time_left <- time_left - execution_time
      } else {
        break
      }
    }
    result[seq_len(length(execution_times))]
  }
}
