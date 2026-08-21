# SPDX-FileCopyrightText: 2026 Stefano Cacciatore
# SPDX-License-Identifier: GPL-2.0-or-later

new_trajectory <- function(input, plotting = FALSE, n = 20, data = NULL,
                           draw = TRUE, knn = 10, trace = NULL, ...) {
  if (plotting) graphics::plot(input, ...)
  x <- input[, 1]
  y <- input[, 2]

  if (is.null(trace)) {
    selected <- graphics::identify(x, y, order = TRUE, plot = FALSE)
    selected <- selected$ind[order(selected$order)]
    spline <- graphics::xspline(
      x[selected], y[selected],
      shape = c(0, rep(-1, length(selected) - 2), 0), draw = FALSE
    )
    new_x <- numeric()
    new_y <- numeric()
    for (i in seq_len(length(spline$x) - 1L)) {
      segment_x <- spline$x[i + 0:1]
      segment_y <- spline$y[i + 0:1]
      prediction_x <- seq(segment_x[1], segment_x[2], length.out = n)
      prediction_y <- stats::predict(
        stats::lm(y ~ x, data = data.frame(x = segment_x, y = segment_y)),
        data.frame(x = prediction_x)
      )
      new_x <- c(new_x, prediction_x)
      new_y <- c(new_y, prediction_y)
    }
    cumulative_distance <- c(
      0,
      cumsum(sqrt(diff(new_x)^2 + diff(new_y)^2))
    )
    sections <- seq(0, utils::tail(cumulative_distance, 1), length.out = n)
    nearest <- Rnanoflann::nn(
      as.matrix(cumulative_distance), as.matrix(sections), 1
    )$indices
    spline$x <- new_x[nearest]
    spline$y <- new_y[nearest]
  } else {
    spline <- trace
  }

  if (draw) {
    graphics::xspline(spline, lwd = 5, border = "white")
    graphics::xspline(spline, lwd = 3, border = "red")
  }
  xy <- cbind(spline$x, spline$y)
  neighbors <- Rnanoflann::nn(input, xy, knn)$indices
  trajectory <- if (is.null(data)) {
    NULL
  } else {
    t(apply(neighbors, 1, function(index) {
      colMeans(as.matrix(data)[index, , drop = FALSE])
    }))
  }
  graphics::points(
    spline$x[1], spline$y[1], col = 2, bg = "#eeeeee", lwd = 2,
    pch = 22, cex = 2
  )
  if (length(spline$x) > 1L) {
    graphics::points(
      spline$x[-1], spline$y[-1], col = 2, bg = "#eeeeee", lwd = 2,
      pch = 21
    )
  }
  list(
    xy = spline, trajectory = trajectory, kk = neighbors,
    settings = list(x = x, y = y, n = n, data = data, knn = knn)
  )
}

volume_rendering <- function(xyz, tissue_segments, selection = NULL,
                             alpha = NULL, colors = NULL,
                             cells = c(20, 20, 20), level = exp(-3)) {
  if (!is.factor(tissue_segments)) stop("tissue_segments is not a factor")
  segments <- levels(tissue_segments)
  if (is.null(colors)) {
    colors <- grDevices::rainbow(length(segments))
  } else if (length(segments) != length(colors)) {
    stop("The number of colors does not match")
  }
  if (is.null(alpha)) {
    alpha <- rep(1, length(segments))
  } else if (length(segments) != length(alpha)) {
    stop("The number of alpha values does not match")
  }
  if (!is.null(selection)) segments <- selection
  selected <- which(levels(tissue_segments) %in% segments)
  colors <- colors[selected]
  alpha <- alpha[selected]
  cells <- pmin(cells, vapply(seq_len(3L), function(j) {
    length(unique(xyz[, j]))
  }, integer(1)))
  visible <- alpha > 0
  segments <- segments[visible]
  alpha <- alpha[visible]
  colors <- colors[visible]

  for (i in seq_along(segments)) {
    selected_segment <- tissue_segments == segments[i]
    density <- misc3d::kde3d(
      xyz[selected_segment, 1], xyz[selected_segment, 2],
      xyz[selected_segment, 3], n = cells
    )
    padded <- array(0, dim = cells + 2)
    padded[2:(cells[1] + 1), 2:(cells[2] + 1), 2:(cells[3] + 1)] <-
      density$d
    dx <- c(
      density$x[1] - density$x[2], density$x,
      density$x[cells[1]] + density$x[2] - density$x[1]
    )
    dy <- c(
      density$y[1] - density$y[2], density$y,
      density$y[cells[2]] + density$y[2] - density$y[1]
    )
    dz <- c(
      density$z[1] - density$z[2], density$z,
      density$z[cells[3]] + density$z[2] - density$z[1]
    )
    misc3d::contour3d(
      padded, level = level, dx, dy, dz, color = colors[i], scale = FALSE,
      engine = "rgl", draw = TRUE, alpha = alpha[i], add = i != 1L
    )
  }
  rgl::rglwidget()
}

plot_slide <- function(xy, slide, labels, col = NULL, nrow = 1,
                       scales = "free", size.dot = 3,
                       size.legend.text = 10, size.legend.title = 10,
                       size.legend.dot = 5, size.strip.text = 10) {
  labels <- as.factor(labels)
  if (is.null(col)) col <- grDevices::rainbow(nlevels(labels))
  frame <- data.frame(
    x = xy[, 1], y = xy[, 2], slide = as.factor(slide), labels = labels
  )
  ggplot2::ggplot(
    frame,
    ggplot2::aes(x = .data$x, y = .data$y, color = .data$labels)
  ) +
    ggplot2::geom_point(size = size.dot) +
    ggplot2::facet_wrap(ggplot2::vars(.data$slide), nrow = nrow, scales = scales) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom", axis.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = size.legend.text),
      legend.title = ggplot2::element_text(size = size.legend.title),
      strip.text = ggplot2::element_text(size = size.strip.text, face = "bold"),
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual("Domain", values = col, drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(
      nrow = 1, override.aes = list(size = size.legend.dot)
    ))
}

read_annotations <- function(address) {
  annotations <- utils::read.csv(address)
  labels <- strsplit(annotations[, 2], ":")
  labels <- unlist(lapply(labels, function(value) value[2]))
  labels <- strsplit(labels, ",")
  labels <- unlist(lapply(labels, function(value) value[1]))
  labels <- gsub('"', "", labels)
  names(labels) <- annotations[, 1]
  occurrence <- stats::ave(
    seq_len(nrow(annotations)), annotations[, 1], FUN = seq_along
  )
  labels <- labels[occurrence == 1]
  substr(labels, 2, nchar(labels))
}
