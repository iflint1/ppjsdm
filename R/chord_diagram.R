#' Plot a chord diagram of the coefficients of a fit object.
#'
#' The coefficients are shown in a chord diagram.
#' Highlighted opaque interactions are statistically significant.
#' The types are ordered according to their mean interaction, or alternatively sorted into classes if the `classes` argument is set.
#'
#' @param fit A fit object obtained by a call to `ppjsdm::gibbsm`.
#' @param coefficient A string representing the coefficient to plot.
#' Choice of `alpha1`, `alpha2`, ..., to show one of the short-range interaction coefficients;
#' `alpha` to show all of the short-range potentials on the same chord diagram;
#' `gamma` to show the medium-range interaction coefficient.
#' @param summ Optional summary corresponding to the fit; if not provided it is obtained by calling `ppjsdm::summary.gibbsm`.
#' @param only_statistically_significant Only show statistically significant coefficients?
#' @param full_names Optional list of full names of types, if for example abbreviations were used when running the fit.
#' @param compute_confidence_intervals Compute the confidence intervals (which is slower) to highlight statistically significant interactions?
#' @param classes If this parameter is supplied, then colours are used to distinguish classes.
#' Should be a named vector/list, with the names corresponding to types, and the value equal to the class name.
#' @param involving Optional vector/list of types. Only coefficients involving these types will be plotted.
#' @param how If the `involving` argument is supplied, should it involve *only* those types, or at least *one* of those types (which is relevant if inter-type interactions are involved).
#' @param include_self Include self-interactions?
#' @param show_legend Show legend(s)?
#' @param cex Cex.
#' @param legend_cex Cex of the legend.
#' @param ninteractions Maximum number of interactions to include.
#' @param track_height Proportion of the chord diagram that the outer rim should occupy.
#' @param show_bottom_left Show the legend in the bottom-left? Only used if `classes` parameter is provided.
#' @param outward_facing_names Should the names of the types be outward facing or along the circle?
#' @param show_grid_ticks Show grid ticks on each of the sectors?
#' @param sort_interactions Should the interactions originating from a given type be sorted?
#' @param circle_margin Margin on the four sides of the circle.
#' @param big_gap Gap size between classes.
#' @param repulsion_attraction_colours Colours to represent repulsion and attraction.
#' @param classes_colours Colours to represent the different classes.
#' @importFrom circlize CELL_META chordDiagramFromDataFrame circos.axis circos.clear circos.initialize circos.link circos.par circos.rect circos.text circos.track circos.trackPlotRegion convert_height get.all.track.index get.cell.meta.data get.current.sector.index get.current.track.index get_most_inside_radius mm_h mm_y rand_color set.current.cell
#' @importFrom graphics legend par
#' @importFrom grDevices adjustcolor col2rgb colorRampPalette rgb
#' @export
#' @md
#' @examples
#' set.seed(1)
#'
#' # Draw a configuration
#' configuration <- ppjsdm::rppp(lambda = c(A = 100, B = 100, C = 100, D = 100))
#'
#' # Fit the data
#' fit <- ppjsdm::gibbsm(configuration)
#'
#' # Chord diagram plot
#' ppjsdm::chord_diagram_plot(fit)
chord_diagram_plot <- function(fit,
                               coefficient = "alpha1",
                               summ,
                               only_statistically_significant = FALSE,
                               full_names = NULL,
                               compute_confidence_intervals = TRUE,
                               classes = NULL,
                               involving = NULL,
                               how = c("only", "one"),
                               include_self = TRUE,
                               show_legend = TRUE,
                               cex = 1,
                               legend_cex = cex,
                               ninteractions = 100,
                               track_height = 0.05,
                               show_bottom_left = TRUE,
                               outward_facing_names = FALSE,
                               show_grid_ticks = TRUE,
                               sort_interactions = TRUE,
                               circle_margin = 0.1,
                               big_gap = 5,
                               repulsion_attraction_colours = c("blue", "red"),
                               classes_colours) {
  # Interpret the how argument
  how <- match.arg(how)

  # Take care of repulsion_attraction_colours argument
  if(!is.character(repulsion_attraction_colours) | length(repulsion_attraction_colours) != 2) {
    stop("repulsion_attraction_colours should be a vector containing two strings, the first one a colour to represent negative interactions, the second one representing positive interactions.")
  }

  # Do we include self-interactions?
  which <- if(include_self) {
    "all"
  } else {
    "between"
  }

  chord_diagram <- make_summary_df(fits = list(fit),
                                   coefficient = coefficient,
                                   summ = summ,
                                   only_statistically_significant = only_statistically_significant,
                                   which = which,
                                   full_names = full_names,
                                   compute_confidence_intervals = compute_confidence_intervals,
                                   classes = classes,
                                   involving = involving,
                                   how = how)

  # Name of the coefficient we are plotting
  identification <- c("gamma", "alpha")
  identification <- identification[identification %in% colnames(chord_diagram)]

  # Remove NAs
  chord_diagram <- chord_diagram[!is.na(df[, identification]), ]

  if(nrow(chord_diagram) == 0) {
    warning("There were no interactions to plot. This could be due to all interactions being non-statistically significant, with option only_statistically_significant active.")
    return()
  }

  if(length(identification) != 1) {
    stop(paste0("Could not identify the intended coefficient by analysing the colnames of the dataframe; identification variable = ", identification))
  }

  # Keep ninteractions largest interactions
  if(ninteractions < nrow(chord_diagram)) {
    chord_diagram <- chord_diagram[order(abs(chord_diagram[, identification]), decreasing = TRUE)[seq_len(ninteractions)], ]
  }

  # Order them from largest interactions to smallest
  chord_diagram <- chord_diagram[order(chord_diagram[, identification], decreasing = TRUE), ]

  # Set colour vector
  col <- ifelse(chord_diagram[, identification] < 0, repulsion_attraction_colours[1], repulsion_attraction_colours[2])
  if(!only_statistically_significant & compute_confidence_intervals) {
    is_significant <- chord_diagram$lo > 0 | chord_diagram$hi < 0
    col <- ifelse(is_significant, adjustcolor(col, alpha.f = 0.8), adjustcolor(col, alpha.f = 0.3))
    link.border <- ifelse(is_significant, "black", "transparent")
  } else {
    col <- adjustcolor(col, alpha.f = 0.8)
    link.border <- "black"
  }

  # Extract type names
  types_names <- unique(c(chord_diagram$from, chord_diagram$to))

  # Compute mean value of interactions for each type
  mean_interaction <- setNames(sapply(types_names, function(ty) {
    mean(chord_diagram[chord_diagram$from == ty | chord_diagram$to == ty, identification])
  }), nm = types_names)

  if(missing(classes)) {
    type_order <- types_names[order(mean_interaction, decreasing = TRUE)]
    grid.col <- setNames(colorRampPalette(c("black", "lightgreen"))(length(types_names)), nm = type_order)
    # gap.after controls the spacing between sectors in the chord diagram
    gap.after <- rep(1, length(type_order))
  } else {
    class_of_types <- setNames(sapply(types_names, function(nm) {
      fr <- chord_diagram$from == nm
      if(sum(fr) > 0) {
        as.character(chord_diagram$class_from[fr][1])
      } else {
        to <- chord_diagram$to == nm
        if(sum(to) > 0) {
          as.character(chord_diagram$class_to[to][1])
        } else {
          stop(paste0("Looking for class of ", nm, " and found nothing informative in the plotting data.frame."))
        }
      }
    }), nm = types_names)
    possible_classes <- sort(unique(class_of_types))

    # Make as many colours as needed, duplicating the supplied ones if needed
    classes_colours <- if(missing(classes_colours)) {
      fixed_colours <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")[seq_len(min(length(possible_classes), 5))]
      setNames(as.list(c(fixed_colours, colorRampPalette(c("black", "lightgreen"))(length(possible_classes) - length(fixed_colours)))),
                       nm = possible_classes)
    } else {
      as.list(classes_colours)
    }

    type_order <- types_names[order(class_of_types, mean_interaction, decreasing = TRUE)]
    grid.col <- setNames(sapply(type_order, function(nm) {
      classes_colours[[class_of_types[nm]]]
    }), nm = type_order)
    # gap.after controls the spacing between sectors in the chord diagram; in this case, we want increased spacing between classes
    gap.after <- rep(1, length(type_order))
    for(cl in unique(grid.col)) {
      last_occ <- which(grid.col == cl)
      last_occ <- last_occ[length(last_occ)]
      gap.after[last_occ] <- big_gap
    }
  }

  # In order to avoid warnings, get rid of superfluous columns
  chord_diagram <- data.frame(from = chord_diagram$from,
                              to = chord_diagram$to,
                              value = chord_diagram[, identification],
                              stringsAsFactors = FALSE)

  if(show_grid_ticks) {
    annotationTrack <- c("grid", "axis")
  } else {
    annotationTrack <- "grid"
  }

  if(sort_interactions) {
    # In this case, the option link.sort from chordDiagramFromDataFrame DOES NOT WORK
    # We therefore use our own patched version (see below) that actually sorts interactions in the way we want
    plot_fun <- .chordDiagramFromDataFrame
  } else {
    plot_fun <- chordDiagramFromDataFrame
  }

  circos.par(gap.after = gap.after, circle.margin = circle_margin)
  plot_fun(chord_diagram,
           reduce = 0, # Default value causes errors when there are tiny sectors
           link.lwd = 1, # Width for link borders
           link.border = link.border, # Solid black lines around links that are highlighted
           col = col,
           self.link = 1, # Self-links as single bumps
           link.zindex = rank(abs(chord_diagram$value)),
           grid.col = grid.col,
           annotationTrack = annotationTrack, # This gets rid of the names on the sectors, these will be printed manually afterwards
           annotationTrackHeight = track_height, # Track size as fraction of total size
           order = type_order,
           preAllocateTracks = 1)

  if(outward_facing_names) {
    facing <- "clockwise"
  } else {
    facing <- "inside"
  }

  if(show_grid_ticks) {
    y_print_text <- CELL_META$ycenter
  } else {
    y_print_text <- CELL_META$ylim[1]
  }

  # Prints the sector names
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(cex = cex,
                x = CELL_META$xcenter,
                y = y_print_text,
                labels = CELL_META$sector.index,
                facing = facing,
                niceFacing = TRUE,
                adj = c(0, 0))
  }, bg.border = NA)

  if(show_legend) {
    if(!missing(classes) & show_bottom_left) {
      legend("bottomleft", inset = 0.01, pch = 15, col = Reduce(c, classes_colours[possible_classes]),
             legend = possible_classes, cex = legend_cex)
    }
    legend("bottomright", inset = 0.01, pch = 15, col = repulsion_attraction_colours,
           legend = c("Repulsion", "Attraction"), cex = legend_cex)
  }

  circos.clear()
}

# Code below is copy/pasted from circlize/R/chordDiagram.R
# MIT @ Zuguang Gu
# The link ordering there does not work in this context.
# We therefore add a few lines after the construction of the data.frame
# used to plot. We rewrite the columns o1, o2, x1 and x2 to order
# links correctly.

.stop_wrap = function(...) {
  x = paste0(...)
  x = paste(strwrap(x), collapse = "\n")
  stop(x, call. = FALSE)
}

.message_wrap = function(...) {
  x = paste0(...)
  x = paste(strwrap(x), collapse = "\n")
  message(x)
}

.warning_wrap = function(...) {
  x = paste0(...)
  x = paste(strwrap(x), collapse = "\n")
  warning(x, call. = FALSE)
}

.parsePreAllocateTracksValue = function(preAllocateTracks) {
  lt = list(ylim = c(0, 1),
            track.height = circos.par("track.height"),
            bg.col = NA,
            bg.border = NA,
            bg.lty = par("lty"),
            bg.lwd = par("lwd"),
            track.margin = circos.par("track.margin"))
  if(length(preAllocateTracks) && is.numeric(preAllocateTracks)) {
    res = vector("list", length = preAllocateTracks)
    for(i in seq_len(preAllocateTracks)) {
      res[[i]] = lt
    }
    return(res)
  } else if(is.list(preAllocateTracks)) {
    # list of list
    if(all(sapply(preAllocateTracks, is.list))) {
      res = vector("list", length = length(preAllocateTracks))
      for(i in seq_along(preAllocateTracks)) {
        lt2 = lt
        for(nm in intersect(names(lt), names(preAllocateTracks[[i]]))) {
          lt2[[nm]] = preAllocateTracks[[i]][[nm]]
        }
        res[[i]] = lt2
      }
      return(res)
    } else {
      lt2 = lt
      for(nm in intersect(names(lt), names(preAllocateTracks))) {
        lt2[[nm]] = preAllocateTracks[[nm]]
      }
      return(list(lt2))
    }
  } else {
    .stop_wrap("Wrong `preAllocateTracks` value.")
  }
}

.convert_h_from_canvas_to_data = function(h,
                                          sector.index = get.current.sector.index(),
                                          track.index = get.current.track.index()) {

  yplot = get.cell.meta.data("yplot", sector.index = sector.index, track.index = track.index)
  cell.ylim = get.cell.meta.data("cell.ylim", sector.index = sector.index, track.index = track.index)
  h/abs(yplot[2] - yplot[1]) * (cell.ylim[2] - cell.ylim[1])
}

.chordDiagramFromDataFrame = function(
    df,
    grid.col = NULL,
    grid.border = NA,
    transparency = 0.5,
    col = NULL,
    order = NULL,
    directional = 0,
    xmax = NULL,
    direction.type = "diffHeight",
    diffHeight = convert_height(2, "mm"),
    link.target.prop = TRUE,
    target.prop.height = mm_h(1),
    reduce = 1e-5,
    self.link = 2,
    preAllocateTracks = NULL,
    annotationTrack = c("name", "grid", "axis"),
    annotationTrackHeight = convert_height(c(3, 2), "mm"),
    link.border = NA,
    link.lwd = par("lwd"),
    link.lty = par("lty"),
    link.auto = TRUE,
    link.sort = "default",
    link.decreasing = TRUE,
    link.arr.length = ifelse(link.arr.type == "big.arrow", 0.02, 0.4),
    link.arr.width = link.arr.length/2,
    link.arr.type = "triangle",
    link.arr.lty = par("lty"),
    link.arr.lwd = par("lwd"),
    link.arr.col = par("col"),
    link.largest.ontop = FALSE,
    link.visible = TRUE,
    link.rank = NULL,
    link.zindex = seq_len(nrow(df)),
    link.overlap = FALSE,
    scale = FALSE,
    group = NULL,
    big.gap = 10,
    small.gap = 1,
    plot = TRUE,
    ...) {

  if(!is.null(link.rank)) {
    .stop_wrap("`link.rank` is not supported and will be removed soon. Please use `link.zindex` instead.")
  }

  if(nrow(df) != 2) {
    if(identical(direction.type, c("diffHeight", "arrows")) || identical(direction.type, c("arrows", "diffHeight"))) {
      direction.type = "diffHeight+arrows"
    }
  }

  # check the format of the data frame
  if(!inherits(df, "data.frame")) {
    .stop_wrap("`df` must be a data frame.")
  }
  df = as.data.frame(df)
  if(ncol(df) < 2) {
    .stop_wrap("`df` should have at least have two columns.")
  }
  if(ncol(df) == 2) {
    df[, 3] = rep(1, nrow(df))
    df[, 4] = df[, 3]
  } else {
    numeric_column = sapply(df, is.numeric); numeric_column[1:2] = FALSE
    if(sum(numeric_column) == 0) {
      df[, 3] = rep(1, nrow(df))
      df[, 4] = df[, 3]
    } else if(sum(numeric_column) == 1) {
      df[, 3] = df[, numeric_column]
      df[, 4] = df[, 3]
    } else if(sum(numeric_column) >= 2) {
      v1 = df[, which(numeric_column)[1]]
      v2 = df[, which(numeric_column)[2]]
      df[, 3] = v1
      df[, 4] = v2
      if(circos.par$message) .message_wrap("There are more than one numeric columns in the data frame. Take the first two numeric columns and draw the link ends with unequal width.\n\nType `circos.par$message = FALSE` to suppress the message.")
    }
  }

  if(scale) {
    for(nm in unique(df[, 1])) df[ df[, 1] == nm, 4] = df[ df[, 1] == nm, 3]/(sum(abs(df[ df[, 1] == nm, 3])) + sum(abs(df[ df[, 2] == nm, 3])))
    for(nm in unique(df[, 2])) df[ df[, 2] == nm, 5] = df[ df[, 2] == nm, 3]/(sum(abs(df[ df[, 2] == nm, 3])) + sum(abs(df[ df[, 1] == nm, 3])))

    df[is.na(df[, 4]), 4] = 0
    df[is.na(df[, 5]), 5] = 0

    df = df[, -3]
  }
  df2 = df[1:4]
  df2[[1]] = as.character(df2[[1]])
  df2[[2]] = as.character(df2[[2]])
  colnames(df2) = c("rn", "cn", "value1", "value2")
  df = df2
  nr = nrow(df)

  if(is.null(link.zindex)) {
    link.zindex = seq_len(nrow(df))
  }

  transparency = ifelse(transparency < 0, 0, ifelse(transparency > 1, 1, transparency))

  cate = union(df[[1]], df[[2]])
  if(!is.null(order)) {

    if(length(grid.col) > 1) {
      if(is.null(names(grid.col))) {
        .warning_wrap("Since you have set `order`, you should better set `grid.col` as a named vector where sector names are the vector names, or else the color will be wrongly assigned.")
      } else {
        if(length(setdiff(cate, names(grid.col))) > 0) {
          .warning_wrap("Since you have set `order`, you should better set `grid.col` as a named vector where sector names are the vector names (should contain all sectors).")
        }
      }
    }

    order = intersect(order, cate)
    if(length(order) != length(cate)) {
      .stop_wrap("`order` should contain names of all sectors.")
    }
    if(is.numeric(order)) {
      if(!setequal(order, seq_along(cate))) {
        .stop_wrap("`order` needs to be integers ranging from 1 to", length(cate))
      }
      cate = cate[order]
    } else {
      if(!setequal(order, cate)) {
        .stop_wrap("`order` should only be picked from sectors.")
      }
      cate = order
    }
  }

  n = length(cate)
  if(is.null(grid.col)) {
    grid.col = rand_color(n)
    names(grid.col) = cate
  } else {
    if(length(grid.col) > 1) {
      if(!is.null(group)) {
        if(is.null(names(grid.col))) {
          .warning_wrap("Since you have set `group`, you should better set `grid.col` as a named vector where sector names are the vector names, or else the color will be wrongly assigned.")
        } else {
          if(length(setdiff(cate, names(grid.col))) > 0) {
            .warning_wrap("Since you have set `group`, you should better set `grid.col` as a named vector where sector names are the vector names (should contain all sectors).")
          }
        }
      }
    }
    if(!is.null(names(grid.col))) {
      unnamed_grid = setdiff(cate, names(grid.col))
      if(length(unnamed_grid) > 0) {
        grid.col = c(grid.col, structure(rand_color(length(unnamed_grid)), names = unnamed_grid))
        # stop("Since your ``grid.col`` is a named vector, all sectors should have corresponding colors.\n")
      }
      grid.col = grid.col[as.vector(cate)]
    } else if(length(grid.col) == 1) {
      grid.col = rep(grid.col, n)
      names(grid.col) = cate
    } else if(length(grid.col) == length(cate)) {
      names(grid.col) = cate
    } else {
      .stop_wrap("Since you set ``grid.col``, the length should be either 1 or number of sectors, or set your ``grid.col`` as vector with names.")
    }
  }

  # colors
  if(is.null(col)) {
    col = grid.col[df[[1]]]
  } else {
    if(is.function(col)) {
      col = col(pmax(df$value1, df$value2))
    } else {
      col = rep(col, nr)[1:nr]
    }
  }

  rgb_mat = t(col2rgb(col, alpha = TRUE))
  if(length(transparency) == nrow(rgb_mat)) {
    if(all(is.na(transparency))) {
      col = rgb(rgb_mat, maxColorValue = 255)
    } else {
      col = rgb(rgb_mat, maxColorValue = 255, alpha = (1-transparency)*255)
    }
  } else {
    if(is.na(transparency)) {
      col = rgb(rgb_mat, maxColorValue = 255, alpha = rgb_mat[, 4])
    } else if(all(rgb_mat[, 4] == 255)) {
      col = rgb(rgb_mat, maxColorValue = 255, alpha = (1-transparency)*255)
    } else {
      col = rgb(rgb_mat, maxColorValue = 255, alpha = rgb_mat[, 4])
    }
  }

  .normalize_to_vector = function(x, link, default) {
    n = nrow(link)
    if(inherits(x, "data.frame")) {
      y = rep(default, n)
      xv = x[[3]]
      names(xv) = paste(x[[1]], x[[2]], sep = "$%^")
      lv = paste(link[[1]], link[[2]], sep = "$%^")
      names(y) = lv
      y[names(xv)] = xv
      y
    } else if(length(x) == 1) {
      rep(x, n)
    } else {
      x
    }
  }

  link.border = .normalize_to_vector(link.border, df[1:2], default = NA)
  link.lwd = .normalize_to_vector(link.lwd, df[1:2], default = 1)
  link.lty = .normalize_to_vector(link.lty, df[1:2], default = 1)
  link.arr.length = .normalize_to_vector(link.arr.length, df[1:2], default = 0.4)
  link.arr.width = .normalize_to_vector(link.arr.width, df[1:2], default = 0.2)
  link.arr.type = .normalize_to_vector(link.arr.type, df[1:2], default = "triangle")
  link.arr.lty = .normalize_to_vector(link.arr.lty, df[1:2], default = 1)
  link.arr.lwd = .normalize_to_vector(link.arr.lwd, df[1:2], default = 1)
  link.arr.col = .normalize_to_vector(link.arr.col, df[1:2], default = NA)
  link.visible = .normalize_to_vector(link.visible, df[1:2], default = NA)
  link.zindex = .normalize_to_vector(link.zindex, df[1:2], default = NA)
  directional = .normalize_to_vector(directional, df[1:2], default = 0)
  direction.type = .normalize_to_vector(direction.type, df[1:2], default = "diffHeight")

  # if(link.auto) {
  # 	od = order(factor(df[, 2], levels = cate),
  # 		       factor(df[, 1], levels = cate))
  # 	df = df[od, , drop = FALSE]
  # 	col = col[od]
  # 	link.border = link.border[od]
  # 	link.lwd = link.lwd[od]
  # 	link.lty = link.lty[od]
  # 	link.arr.length = link.arr.length[od]
  # 	link.arr.width = link.arr.width[od]
  # 	link.arr.type = link.arr.type[od]
  # 	link.arr.lty = link.arr.lty[od]
  # 	link.arr.lwd = link.arr.lwd[od]
  # 	link.arr.col = link.arr.col[od]
  # 	link.visible = link.visible[od]
  # 	link.zindex = link.zindex[od]
  # 	directional = directional[od]
  # 	direction.type = direction.type[od]
  # }

  #### reduce the data frame
  onr = nrow(df)
  onn = union(df[, 1], df[, 2])
  while(1) {
    xsum = structure(rep(0, length(cate)), names = cate)
    for(i in seq_len(nr)) {
      if(df$rn[i] == df$cn[i]) {
        xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value1[i])
        if(self.link == 2) {
          xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value2[i])  # <<- self-link!!!!!
        }
      } else {
        xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value1[i])
        xsum[df$cn[i]] = xsum[df$cn[i]] + abs(df$value2[i])
      }
    }

    keep = names(xsum)[xsum / sum(xsum) >= reduce]
    l = df$rn %in% keep & df$cn %in% keep

    cate = intersect(cate, keep)
    df = df[l, , drop = FALSE]
    grid.col = grid.col[intersect(names(grid.col), keep)]
    col = col[l]
    link.border = link.border[l]
    link.lwd = link.lwd[l]
    link.lty = link.lty[l]
    link.arr.length = link.arr.length[l]
    link.arr.width = link.arr.width[l]
    link.arr.type = link.arr.type[l]
    link.arr.lwd = link.arr.lwd[l]
    link.arr.lty = link.arr.lty[l]
    link.arr.col = link.arr.col[l]
    link.visible = link.visible[l]
    link.zindex = link.zindex[l]
    directional = directional[l]
    direction.type = direction.type[l]

    nr = nrow(df)
    reduce = 1e-10
    if(nr == onr) break
    onr = nr
  }

  nn = union(df[, 1], df[, 2])
  if(length(nn) < length(onn)) {
    gap.degree = circos.par$gap.degree
    if(length(gap.degree) > 1) {
      if(is.null(names(gap.degree))) {
        if(length(nn) != length(gap.degree)) {
          .stop_wrap("`reduce` argument causes reduction of sectors. You can either set `reduce = 0`, or adjust the `gap.degree`/`gap.after` in `circos.par()` to remove tiny sectors, or set `gap.degree`/`gap.after` as a named vector where sector names are the vector names (tiny sectors can be included).")
        }
      } else {
        if(length(setdiff(nn, names(gap.degree))) == 0) {

        } else {
          if(length(nn) != length(gap.degree)) {
            .stop_wrap("`reduce` argument causes reduction of sectors. You can either set `reduce = 0`, or adjust the `gap.degree`/`gap.after` in `circos.par()` to remove tiny sectors, or set `gap.degree`/`gap.after` as a named vector where sector names are the vector names (tiny sectors can be included).")
          }
        }
      }
    }
  }

  # re-calcualte xsum
  xsum = structure(rep(0, length(cate)), names = cate)
  for(i in seq_len(nr)) {
    if(df$rn[i] == df$cn[i]) {
      xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value1[i])
      if(self.link == 2) {
        xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value2[i])  # <<- self-link!!!!!
      }
    } else {
      xsum[df$rn[i]] = xsum[df$rn[i]] + abs(df$value1[i])
      xsum[df$cn[i]] = xsum[df$cn[i]] + abs(df$value2[i])
    }
  }

  ### group_by ###
  gap.after = NULL
  if(!is.null(group)) {
    # validate `group`
    if(is.null(names(group))) {
      .stop_wrap("`group` should be named vector where names are the sector names and values are the group labels.")
    }
    if(any(duplicated(names(group)))) {
      .stop_wrap("Names in `group` should not be duplicated.")
    }
    if(length(setdiff(cate, names(group))) > 0) {
      .stop_wrap("Names in `group` should cover all sector names.")
    }

    group = group[intersect(names(group), cate)]

    tg = table(group)
    group_lt = split(names(group), group)

    sn_by_group = unlist(group_lt)

    cate = intersect(sn_by_group, cate)
    xsum = xsum[cate]

    gap.after = c(unlist(lapply(tg, function(x) c(rep(small.gap, x-1), big.gap))))
    names(gap.after) = sn_by_group
  }

  # add additinal columns in df
  df$o1 = rep(0, nr)  # order of the link root in the sector
  df$o2 = rep(0, nr)  # order of the other link root in the sector
  df$x1 = rep(0, nr)  # end position of the link root in the sector
  df$x2 = rep(0, nr)  # end position on the other link root in the sector

  ######## sort links on every sector
  # row first
  .order = function(x, sort = FALSE, reverse = FALSE) {
    if(sort) {
      od = order(x)
    } else {
      od = seq_along(x)
    }
    if(reverse) {
      od = rev(od)
    }
    return(od)
  }

  if(length(link.sort) == 1) link.sort = rep(link.sort, 2)
  if(length(link.decreasing) == 1) link.decreasing = rep(link.decreasing, 2)

  if(identical(link.sort[1], "default")) {
    od = tapply(seq_len(nrow(df)), df$rn, function(ind) {
      rn = df[ind[1], "rn"]
      cn = df[ind, "cn"]
      fa = c(cate, cate)
      i = which(fa == rn)[1]
      fa = fa[seq(i, i + length(cate) - 1)]
      order(factor(cn, levels = fa), decreasing = TRUE)
    })
  } else if(identical(link.sort[1], "asis")) {
    od = tapply(abs(df$value1), df$rn, .order, FALSE, FALSE)
  } else if(identical(link.sort[1], FALSE)) {
    od = tapply(abs(df$value1), df$rn, .order, FALSE, link.decreasing[1])
  } else {
    # position of root 1
    od = tapply(abs(df$value1), df$rn, .order, link.sort[1], link.decreasing[1])
  }

  for(nm in names(od)) {  # for each sector
    l = df$rn == nm # rows in df that correspond to current sector
    df$o1[l] = od[[nm]] # adjust rows according to the order in current sector
    df$x1[l][od[[nm]]] = cumsum(abs(df$value1[l])[od[[nm]]]) # position

    l2 = df$rn == nm & df$cn == nm
    if(sum(l2)) { # there is a self link
      if(self.link == 1) {
        df$x2[l2] = df$x1[l2]+abs(df$value1[l2])*0.000001
      }
    }
  }
  max_o1 = sapply(od, max)
  sum_1 = tapply(abs(df$value1), df$rn, sum)
  # position of root 2
  if(identical(link.sort[2], "default")) {
    od = tapply(seq_len(nrow(df)), df$cn, function(ind) {
      cn = df[ind[1], "cn"]
      rn = df[ind, "rn"]
      fa = c(cate, cate)
      i = which(fa == cn)[1]
      fa = fa[seq(i, i + length(cate) - 1)]
      # for the ordering of duplicated rows
      v = numeric(length(rn))
      for(x in unique(rn)) {
        l = rn == x
        v[l] = seq_len(sum(l))
      }
      order(factor(rn, levels = fa), v, decreasing = TRUE)
    })
  } else if(identical(link.sort[2], "asis")) {
    od = tapply(abs(df$value2), df$cn, .order, FALSE, FALSE)
  } else if(identical(link.sort[2], FALSE)) {
    od = tapply(abs(df$value2), df$cn, .order, FALSE, link.decreasing[2])
  } else {
    od = tapply(abs(df$value2), df$cn, .order, link.sort[2], link.decreasing[2])
  }

  for(nm in names(od)) {
    if(!is.na(max_o1[nm])) { # if cn already in rn
      l = df$cn == nm
      if(self.link == 1) {
        l2 = ! df$rn[l] == nm # self link
        # od[[nm]][l2] = order(od[[nm]][l2])
        if(sum(l2)) {
          od[[nm]] = order(od[[nm]][l2])
        }
      } else {
        l2 = rep(TRUE, sum(l))
      }
      df$o2[l][l2] = od[[nm]] + max_o1[nm]
      df$x2[l][l2][ od[[nm]] ] = cumsum(abs(df$value2[l][l2])[ od[[nm]] ]) + sum_1[nm]
    } else {
      l = df$cn == nm
      df$o2[l] = od[[nm]]
      df$x2[l][od[[nm]]] = cumsum(abs(df$value2[l])[od[[nm]]])
    }
  }
  if(self.link == 1) {
    l = df$rn == df$cn
    df$x1[l] = pmin(df$x1[l], df$x2[l])
    df$x2[l] = pmin(df$x1[l], df$x2[l])
  }

  ################################################################
  # Part changed is between these comments
  ################################################################

  type_list <- unique(union(df$rn, df$cn))
  for(ty in type_list) {
    contains_ty <- df$cn == ty | df$rn == ty
    ord <- order(df$value1[contains_ty], decreasing = link.decreasing[1])
    k <- 1
    val <- 0.
    for(j in seq_len(nrow(df))) {
      if(df$rn[j] == ty) {
        df$o1[j] <- ord[k]
        val <- val + abs(df$value1[j])
        df$x1[j] <- val
        k <- k + 1
      } else if(df$cn[j] == ty) {
        val <- val + abs(df$value1[j])
        k <- k + 1
      }
    }
    ord <- order(df$value2[contains_ty], decreasing = link.decreasing[2])
    k <- 1
    val <- 0.
    for(j in seq_len(nrow(df))) {
      if(df$cn[j] == ty) {
        df$o2[j] <- ord[k]
        val <- val + abs(df$value2[j])
        df$x2[j] <- val
        k <- k + 1
      } else if(df$rn[j] == ty) {
        val <- val + abs(df$value2[j])
        k <- k + 1
      }
    }
  }

  ################################################################
  # Part changed is between these comments
  ################################################################

  if(!is.null(xmax)) {
    overlap = intersect(names(xmax), names(xsum))
    xmax = xmax[overlap]
    xmax = xmax[xmax > xsum[overlap]]
    if(length(xmax)) {
      xsum[names(xmax)] = xmax
    }
  }

  if(link.overlap) {

    if(!directional) {
      warning("`link.overlap` should be used with directional = 1 or -1.")
    }

    x1_sum = tapply(df$x1, df$rn, max)[names(xsum)]
    x1_sum[is.na(x1_sum)] = 0; names(x1_sum) = names(xsum)
    x2_sum = tapply(df$x2, df$cn, max)[names(xsum)] - x1_sum
    x2_sum[is.na(x2_sum)] = 0; names(x2_sum) = names(xsum)
    df$x2 = df$x2 - x1_sum[df$cn]
    xsum = pmax(x1_sum, x2_sum)
  }

  if(!plot) return(df)

  o.cell.padding = circos.par("cell.padding")
  circos.par(cell.padding = c(0, 0, 0, 0))
  o.start.degree = circos.par("start.degree")
  o.gap.after = circos.par("gap.after")
  o.points.overflow.warning = circos.par("points.overflow.warning")
  circos.par("points.overflow.warning" = FALSE)
  ### group_by ###
  if(!is.null(gap.after)) {
    circos.par(gap.after = gap.after)
  }

  if((identical(circos.par("gap.after"), 1))) {  # gap.after is not set in circos.par()
    if(length(intersect(df[, 1], df[, 2])) == 0) {
      ind1 = which(cate %in% df[, 1])
      ind2 = which(cate %in% df[, 2])

      if(max(ind1) < min(ind2) || max(ind2) < min(ind1)) {
        n_df1 = length(unique(df[, 1]))
        n_df2 = length(unique(df[, 2]))
        s1 = sum(abs(df[, 3]))
        s2 = sum(abs(df[, 4]))
        if(cate[1] %in% df[ ,2]) {
          foo = n_df1
          n_df1 = n_df2
          n_df2 = foo

          foo = s1
          s1 = s2
          s2 = foo
        }
        n_sector = n_df1 + n_df2
        d1 =  (360 - small.gap*(n_sector - 2) - big.gap*2) * (s1/(s1 + s2)) + small.gap*(n_df1-1)
        if(circos.par$start.degree == 1) circos.par$start.degree = 0
        if(circos.par$clock.wise) {
          start_degree = circos.par$start.degree - (180 - d1)/2
        } else {
          start_degree = circos.par$start.degree + (180 - d1)/2
        }
        gap.after = c(rep(small.gap, n_df1 - 1), big.gap, rep(small.gap, n_df2 - 1), big.gap)
        suppressWarnings(circos.par(start.degree = start_degree, gap.after = gap.after))
      }
    } else {
      # warning("The two sets of sectors overlap, ignore `gap.degree`.")
    }
  } else {
    # warning("You have changed the default value of circos.par('gap.degree') or circos.par('gap.after').\nIgnore `gap.degree` argument.")
  }
  circos.initialize(factors = factor(cate, levels = cate), xlim = cbind(rep(0, length(xsum)), xsum))

  # pre allocate track
  if(!is.null(preAllocateTracks)) {
    pa = .parsePreAllocateTracksValue(preAllocateTracks)
    for(i in seq_along(pa)) {
      va = pa[[i]]
      circos.trackPlotRegion(ylim = va$ylim, track.height = va$track.height,
                             bg.col = va$bg.col, bg.border = va$bg.border, bg.lty = va$bg.lty, bg.lwd = va$bg.lwd,
                             track.margin = va$track.margin)
    }
  }
  if("name" %in% annotationTrack) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA,
                           panel.fun = function(x, y) {
                             xlim = get.cell.meta.data("xlim")
                             current.sector.index = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             circos.text(mean(xlim), 0.9, labels = current.sector.index, cex = 0.8,
                                         facing = "inside", niceFacing = TRUE, adj = c(0.5, 0))
                           }, track.height = annotationTrackHeight[which(annotationTrack %in% "name")])
  }
  if("grid" %in% annotationTrack) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA,
                           panel.fun = function(x, y) {
                             xlim = get.cell.meta.data("xlim")
                             current.sector.index = get.cell.meta.data("sector.index")
                             if(is.null(grid.border)) {
                               border.col = grid.col[current.sector.index]
                             } else {
                               border.col = grid.border
                             }
                             circos.rect(xlim[1], 0, xlim[2], 1, col = grid.col[current.sector.index], border = border.col)
                             if("axis" %in% annotationTrack) {
                               if(scale) {
                                 circos.axis("top", labels = function(x) {paste0(round(x*100), "%")}, labels.cex = 0.5)
                               } else {
                                 circos.axis("top", labels.cex = 0.5)
                               }
                             }
                           }, track.height = annotationTrackHeight[which(annotationTrack %in% "grid")])
  }

  rou = get_most_inside_radius()
  rou1 = numeric(nr)
  rou2 = numeric(nr)
  for(i in seq_len(nr)) {
    if(directional[i]) {
      if(grepl("diffHeight", direction.type[i])) {
        if(directional[i] == 1) {
          if(diffHeight > 0) {
            rou1[i] = rou - diffHeight
            rou2[i] = rou
          } else {
            rou1[i] = rou
            rou2[i] = rou + diffHeight
          }
        } else if(directional[i] == -1) {
          if(diffHeight > 0) {
            rou1[i] = rou
            rou2[i] = rou - diffHeight
          } else {
            rou1[i] = rou + diffHeight
            rou2[i] = rou
          }
        } else  if(directional[i] == 2) {
          if(diffHeight > 0) {
            rou1[i] = rou - diffHeight
            rou2[i] = rou - diffHeight
          } else {
            rou1[i] = rou + diffHeight
            rou2[i] = rou + diffHeight
          }
        }
      } else {
        rou1[i] = rou
        rou2[i] = rou
      }
    } else {
      rou1[i] = rou
      rou2[i] = rou
    }
  }

  if(link.largest.ontop) {
    link_order = order(pmax(abs(df$value1), abs(df$value2)), decreasing = FALSE)
  } else {
    link_order = order(link.zindex)
  }

  for(k in link_order) {
    if(abs(df$value1[k])/sum(abs(df$value1)) < 1e-6 && abs(df$value2[k])/sum(abs(df$value2)) < 1e-6) next
    if(link.visible[k] && col[k] != "#FFFFFF00") {
      if(setequal(direction.type, c("diffHeight"))) {
        circos.link(df$rn[k], c(df$x1[k] - abs(df$value1[k]), df$x1[k]),
                    df$cn[k], c(df$x2[k] - abs(df$value2[k]), df$x2[k]),
                    directional = 0, col = col[k], rou1 = rou1[k], rou2 = rou2[k],
                    border = link.border[k], lwd = link.lwd[k], lty = link.lty[k],
                    ...)
      } else if(any(grepl("arrows", direction.type))) {
        circos.link(df$rn[k], c(df$x1[k] - abs(df$value1[k]), df$x1[k]),
                    df$cn[k], c(df$x2[k] - abs(df$value2[k]), df$x2[k]),
                    directional = directional[k], col = col[k], rou1 = rou1[k], rou2 = rou2[k],
                    border = link.border[k],
                    lwd = link.lwd[k], lty = link.lty[k],
                    arr.length = link.arr.length[k], arr.width = link.arr.width[k],
                    arr.type = link.arr.type[k], arr.col = link.arr.col[k],
                    arr.lty = link.arr.lty[k], arr.lwd = link.arr.lwd[k],
                    ...)
      }
    }
  }

  df$col = col

  if(link.target.prop && all(grepl("diffHeight", direction.type)) && min(diffHeight) > 0) {
    if(all(directional %in% c(1, 2))) {
      last.track.index = rev(get.all.track.index())[1]
      for(i in seq_len(nrow(df))) {
        if(abs(df$value1[i]) > 0 && link.visible[i] && col[i] != "#FFFFFF00") {
          set.current.cell(sector.index = df$rn[i], track.index = last.track.index)
          bar_h = .convert_h_from_canvas_to_data(min(target.prop.height, min(diffHeight)))
          rect_y = .convert_h_from_canvas_to_data(get.cell.meta.data("track.margin", track.index = last.track.index)[1] + circos.par("track.margin")[2])
          circos.rect(df[i, "x1"], -rect_y,
                      df[i, "x1"] - abs(df[i, "value1"]), -rect_y - bar_h,
                      col = grid.col[df$cn[i]], border = NA)
        }
      }
    }
    if(all(directional %in% c(-1, 2))) {
      last.track.index = rev(get.all.track.index())[1]
      for(i in seq_len(nrow(df))) {
        if(abs(df$value1[i]) > 0 && link.visible[i] && col[i] != "#FFFFFF00") {
          set.current.cell(sector.index = df$cn[i], track.index = last.track.index)
          bar_h = .convert_h_from_canvas_to_data(min(target.prop.height, min(diffHeight)))
          circos.rect(df[i, "x2"], -mm_y(1),
                      df[i, "x2"] - abs(df[i, "value2"]), -mm_y(1) - bar_h,
                      col = grid.col[df$rn[i]], border = NA)
        }
      }
    }
  }

  suppressWarnings(circos.par("cell.padding" = o.cell.padding, "start.degree" = o.start.degree,
                              "gap.after" = o.gap.after, "points.overflow.warning" = o.points.overflow.warning))
  return(invisible(df))
}

