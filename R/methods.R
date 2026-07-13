#' @export
coef.pdmi_frailty <- function(object, ...) {
  pdmi_result_omega(object)
}

#' @export
vcov.pdmi_frailty <- function(object, ...) {
  pdmi_result_vcov(object)
}

#' Print a PDMI frailty fit
#'
#' @param x A `pdmi_frailty` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.pdmi_frailty <- function(x, ...) {
  cat("pdmi_frailty fit\n")
  cat("  missing type: ", x$missing_type, "\n", sep = "")
  cat("  prior version: ", x$prior_version, "\n", sep = "")
  cat("  imputations: ", x$M, "; burn-in: ", x$B,
      "; retained updates: ", x$numit, "\n", sep = "")
  if (isTRUE(x$convergence)) {
    cat("  convergence: TRUE\n")
  } else {
    cat("  convergence: FALSE\n")
    cat("  failure reason: ", x$failure_reason %||% "unknown", "\n", sep = "")
  }
  if (x$report %in% c("both", "parameters")) {
    tab <- pdmi_parameter_table(x)
    print(tab, row.names = FALSE, digits = 4)
  }
  if (x$report %in% c("both", "penetrance") &&
      !is.null(x$result$penetrance) && nrow(x$result$penetrance)) {
    cat("  penetrance rows: ", nrow(x$result$penetrance), "\n", sep = "")
  }
  invisible(x)
}

#' Summarize pooled PDMI parameter estimates
#'
#' @param object A `pdmi_frailty` object.
#' @param conf.level Confidence level for Wald intervals.
#' @param ... Unused.
#'
#' @return A `summary.pdmi_frailty` object.
#' @export
summary.pdmi_frailty <- function(object,
                                 conf.level = object$conf.level,
                                 ...) {
  out <- list(
    call = object$call,
    convergence = object$convergence,
    failure_reason = object$failure_reason,
    missing_type = object$missing_type,
    prior_version = object$prior_version,
    coefficients = pdmi_parameter_table(object, conf.level = conf.level),
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.pdmi_frailty"
  out
}

#' @export
print.summary.pdmi_frailty <- function(x, ...) {
  cat("pdmi_frailty parameter summary\n")
  cat("  missing type: ", x$missing_type,
      "; prior version: ", x$prior_version, "\n", sep = "")
  if (!isTRUE(x$convergence)) {
    cat("  convergence: FALSE\n")
    cat("  failure reason: ", x$failure_reason %||% "unknown", "\n", sep = "")
  }
  print(x$coefficients, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Return the derived congenial imputation model
#'
#' @description
#' Prints and returns the kernel components that `pdmi_frailty()` derived from
#' the analysis formula and `prior`. This is the object-level answer to whether
#' the congenial imputation model was generated automatically.
#'
#' @param x A fitted object.
#' @param ... Unused.
#'
#' @return A `pdmi_imputation_model` object.
#' @export
imputation_model <- function(x, ...) {
  UseMethod("imputation_model")
}

#' @export
imputation_model.pdmi_frailty <- function(x, ...) {
  print(x$model)
  invisible(x$model)
}

#' @export
print.pdmi_imputation_model <- function(x, ...) {
  cat("Derived frailtypack-congenial PDMI kernel\n")
  cat("  target: ", x$technical_note_target, "\n", sep = "")
  cat("  automatic: ", x$automatic, "\n", sep = "")
  cat("  missing type: ", x$missing_type,
      "; prior version: ", x$prior_version, "\n", sep = "")
  cat("  analysis formula: ", x$analysis_formula, "\n", sep = "")
  cat("  kernel: ", x$kernel, "\n", sep = "")
  cat("  disease: ", x$disease_component, "\n", sep = "")
  cat("  ascertainment: ", x$ascertainment_component, "\n", sep = "")
  cat("  support: ", x$support_component, "\n", sep = "")
  cat("  sampler:\n")
  cat("    theta: ", x$sampler$disease_parameter_draw, "\n", sep = "")
  cat("    frailty: ", x$sampler$frailty_draw, "\n", sep = "")
  cat("    continuous: ", x$sampler$continuous_missing_draw, "\n", sep = "")
  cat("    binary: ", x$sampler$binary_missing_draw, "\n", sep = "")
  cat("    joint: ", x$sampler$joint_extra_term, "\n", sep = "")
  invisible(x)
}

#' Summarize penetrance estimates
#'
#' @param x A `pdmi_frailty` object.
#' @param ... Unused.
#'
#' @return A `pdmi_penetrance_summary` data frame.
#' @export
pen_summary <- function(x, ...) {
  UseMethod("pen_summary")
}

#' @rdname pen_summary
#' @param conf.level Confidence level for intervals.
#' @param ci If `TRUE`, include `lower` and `upper`.
#' @export
pen_summary.pdmi_frailty <- function(x,
                                     conf.level = x$conf.level,
                                     ci = x$pen_ci,
                                     ...) {
  pen <- x$result$penetrance
  if (is.null(pen) || !nrow(pen)) {
    stop("No penetrance estimates are available in this fitted object.",
         call. = FALSE)
  }
  out <- as.data.frame(pen)
  if (isTRUE(ci)) {
    z <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$lower <- pmax(0, out$estimate - z * out$se)
    out$upper <- pmin(1, out$estimate + z * out$se)
  }
  attr(out, "conf.level") <- conf.level
  class(out) <- c("pdmi_penetrance_summary", "data.frame")
  out
}

#' @export
print.pdmi_penetrance_summary <- function(x, ...) {
  cat("pdmi_frailty penetrance summary")
  cl <- attr(x, "conf.level")
  if (!is.null(cl)) cat(" (", 100 * cl, "% CI)", sep = "")
  cat("\n")
  print.data.frame(x, row.names = FALSE, digits = 4)
  invisible(x)
}
