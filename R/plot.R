#' Plot penetrance curves
#'
#' @param x A `pdmi_frailty` object.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
pen_plot <- function(x, ...) {
  UseMethod("pen_plot")
}

#' @rdname pen_plot
#' @param ci If `TRUE`, draw confidence interval ribbons when available.
#' @param conf.level Confidence level used when intervals need to be added.
#' @export
pen_plot.pdmi_frailty <- function(x,
                                  ci = x$pen_ci,
                                  conf.level = x$conf.level,
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required for pen_plot(). Install ggplot2 or use pen_summary().",
         call. = FALSE)
  }
  tab <- pen_summary(x, conf.level = conf.level, ci = ci)
  tab$profile <- paste0("PRS ", tab$prs, ", carrier ", tab$gene)
  p <- ggplot2::ggplot(
    tab,
    ggplot2::aes(x = age, y = estimate, color = profile, group = profile)
  )
  if (isTRUE(ci) && all(c("lower", "upper") %in% names(tab))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = profile),
      alpha = 0.16,
      color = NA
    )
  }
  p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::facet_grid(gene ~ prs, labeller = ggplot2::label_both) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Age",
      y = "Penetrance",
      color = "Profile",
      fill = "Profile"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
