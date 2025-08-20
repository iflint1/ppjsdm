#' Heatmap of coefficient values of a fit object
#'
#' @param fit A fit object obtained by a call to 'ppjsdm::gibbsm'
#' @param summ Provide a summary object of the fit, only needed if show_sig = TRUE
#' @param coefficient A string representing the coefficient to plot, either alpha or beta
#' @param show_values Should the cells show the coefficient value?
#' @param show_sig Should the significance of each cell be shown by a point in the cell?
#' @param limits Optional vector that gives what coefficient values should be plotted
#' @param colours Colours that the gradient of the heatmap will be filled by
#' @param scale_by For visualisation - what should the interval of legend scale increase by?
#' @param include_y Optional vector of types that should be included on the y-axis
#' @param include_x Optional vector of types that should be included on the x-axis
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit
#' @importFrom ggplot2 aes geom_tile scale_fill_gradientn labs theme geom_text geom_point scale_size_manual waiver
#' @importFrom scales squish
#' @importFrom rlang splice-operator
#' @importFrom dplyr filter mutate recode
#' @returns
#' @export
#'
#' @examples
#' #' set.seed(1)
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
                     summ = NULL,
                     coefficient = c("alpha", "beta"),
                     show_values = TRUE,
                     show_sig = FALSE,
                     limits = NULL,
                     colours = c("#6c63ac", "#61AEE6","#80D6BE","gray95", "#F2CE61","#F5A54F","#F55E2C"),
                     scale_by = 0.5,
                     include_y = NULL,
                     include_x = NULL,
                     full_names = NULL)
{

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

if(!is.null(summ)){
#insert low CI value
df$lo <- sapply(seq_len(nrow(df)), function(i) { # Get the lower-bound of the CIs
  val <- summ$lo$alpha[[1]][df$from[i], df$to[i]]
  if(length(val) == 0) {
    val <- summ$lo$alpha[[1]][df$to[i], df$from[i]]
  }
  val
})

#insert high CI value
df$hi <- sapply(seq_len(nrow(df)), function(i) { # Get the lower-bound of the CIs
  val <- summ$hi$alpha[[1]][df$from[i], df$to[i]]
  if(length(val) == 0) {
    val <- summ$hi$alpha[[1]][df$to[i], df$from[i]]
  }
  val
})

df$sig <- sapply(seq_len(nrow(df)), function(i) {
  sig <- ifelse(df$lo[i] > 0 | df$hi[i] < 0, "1", "0") #make significance column, 1 = significant, 0 = not significant
sig
  })

}


if(!is.null(include_x)){

  df <- df %>% filter(from %in% include_x)

}

if(!is.null(include_y)){

  df <- df %>% filter(to %in% include_y)

}

if(!is.null(full_names)){
  df <- df %>% mutate(from = recode(from, !!!names))

  df <- df %>% mutate(to = recode(to, !!!names))
}


if(!is.null(limits)){
breaks <- seq(limits[1], limits[2], by = scale_by)

custom_labels <- c(paste0(">", limits[1]), breaks[2:(length(breaks) - 1)], paste0(">",limits[2]))
} else{
  breaks <- waiver()
  custom_labels <- waiver()
}


p <- ggplot(df, aes(x = from,
                    y = to,
                    fill = value)) +

     geom_tile(colour = "white",
               height = 1.0,
               alpha = 0.75 ) +

     scale_fill_gradientn(colours = colours,
                          limits = limits,
                          breaks = breaks,
                          labels = custom_labels,
                          name = "Coefficient",
                          oob = scales::squish) +

     labs(x = "",
          y = "") +

    theme(panel.grid.major = element_blank(),

          panel.grid.minor = element_blank(),

          panel.background = element_blank(),

         legend.position = "right",

         legend.title=element_text(size=15),

          legend.text=element_text(size=12),

         axis.title.x = element_text(size=25),

          axis.title.y = element_text(size=25),

          axis.text.x = element_text(size=15, angle = 45, hjust = 1, face = "bold"),

        axis.text.y = element_text(size=15, face = "bold"),

          strip.text.y = element_text(size = 15),

          strip.text.x = element_text(size = 15))




if(show_values == TRUE){
  p <- p + geom_text(aes(label = round(value, 2)), color = "black", size = 4.25)
  }

if(show_sig == TRUE) {

df$sig <- factor(df$sig)

  p <- p +
    geom_point(aes(size = as.factor(sig))) +
    scale_size_manual(values = c(NA, 2),
                      name = " ",
                      breaks = c("", "1"),
                      labels = c("", "Significant"))

  }

 print(p)

 }

