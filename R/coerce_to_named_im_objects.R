add_names <- function(str, covariates) {
  if(is.null(names(covariates))) {
    no_name <- rep(TRUE, length(covariates))
  } else {
    no_name <- names(covariates) == ""
  }
  names(covariates)[no_name] <- sprintf(paste0(str, "%d"), seq_len(length(which(no_name))))
  covariates
}

coerce_to_named_im_objects <- function(lst, str, window) {
  lst <- lapply(as.list(lst), function(element) as.im(element, W = as.owin(window)))
  add_names(str, lst)
}
