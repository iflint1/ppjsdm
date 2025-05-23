#' Heatmap of coefficient values of a fit object
#'
#' @param fit A fit object obtained by a call to 'ppjsdm::gibbsm'
#' @param coefficient A string representing the coefficient to plot
#' @param show_values Should the cells show the coefficient value?
#' @param limits Optional vector that gives what coefficient values should be plotted
#' @param include_y Optional vector of types that should be included on the y-axis
#' @param include_x Optional vector of types that should be included on the x-axis
#' @param names Optional list of full names of types, if for example abbreviations were used when running the fit
#' @importFrom dplyr filter mutate recode
#' @importFrom ggplot2 aes geom_tile scale_fill_gradient2 labs theme geom_text waiver
#' @importFrom rlang !!!
#' @importFrom scales squish
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' # Construct a configuration
#'
#' configuration <- ppjsdm::rppp(lambda = c(A = 500, B = 500))
#'
#' # Fit the model
#'
#' fit <- ppjsdm::gibbsm(configuration, covariates = list(x = function(x, y) x))
#'
#' # Plot the heatmap of coefficients
#'
#' ppjsdm::heat_map(fit, coefficient = "x", show_values = FALSE)
#'
heat_map <- function(fit,
                     coefficient = c("alpha", "beta"),
                     show_values = TRUE,
                     limits = NULL,
                     include_y = NULL,
                     include_x = NULL,
                     names = NULL) {

  if(coefficient == "alpha") {
    estimates <- fit$coefficients$alpha[[1]]
  } else {
    estimates <- fit$coefficients$beta}


  unique_names <- colnames(estimates)
  df <- as.data.frame(expand.grid(from = rownames(estimates), to = colnames(estimates)))
  #insert coefficient values
  df$value <- sapply(seq_len(nrow(df)), function(i) {
    val <- estimates[df$from[i], df$to[i]]
    if(length(val) == 0) {
      val <- estimates[df$to[i], df$from[i]]
    }
    val
  })

  if(!is.null(include_x)){

    df <- df %>% filter(from %in% include_x)

  }

  if(!is.null(include_y)){

    df <- df %>% filter(to %in% include_y)

  }

  if(!is.null(names)){
    df <- df %>% mutate(from = recode(from, !!!names))

    df <- df %>% mutate(to = recode(to, !!!names))
  }


  if(!is.null(limits)){
    breaks <- seq(limits[1], limits[2], by = 1)

    custom_labels <- c(paste0(">", limits[1]), breaks[2:(length(breaks) - 1)], paste0(">", limits[2]))
  } else{
    breaks <- waiver()
    custom_labels <- waiver()
  }


  p <- ggplot(df, aes(x = from, y = to, fill = value)) +
    geom_tile(colour = "white", height = 1.0, alpha = 0.75 ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = limits, breaks = breaks, labels = custom_labels, name = "Coefficient", oob = squish) +
    labs(x = "", y = "") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "right",
          legend.title=element_text(size = 15),
          legend.text=element_text(size = 12),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 15, angle = 45, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(size = 15),
          strip.text.x = element_text(size = 15))




  if(show_values == TRUE){
    p <- p + geom_text(aes(label = round(value, 2)), color = "black", size = 4.25)
  }


  print(p)

}

