#' Control settings for `pdmi_frailty()`
#'
#' @description
#' Builds a list of numerical and sampler controls used by the V2.2 exact-slice
#' PDMI engine. Arguments not listed here can be passed through `...` when they
#' correspond to a V2.2 configuration field.
#'
#' @param agemin Minimum analysis age used in the Weibull baseline and
#'   ascertainment denominator.
#' @param frailtypack_maxit Maximum number of frailtypack iterations for each
#'   completed-data model fit.
#' @param gh_order Gauss-Hermite quadrature order. V2.2 supports 10 or 20.
#' @param theta_slice_widths Named slice widths for `log.rho`, `log.lambda`,
#'   `beta_b`, `beta_c`, and `sigma_u2`.
#' @param theta_slice_sweeps Number of parameter slice sweeps per PDMI update.
#' @param theta_slice_m Maximum stepping-out steps for disease-parameter slice
#'   updates.
#' @param continuous_slice_step Slice step size retained for compatibility with
#'   earlier V2 code.
#' @param continuous_slice_m Maximum continuous-sampler stepping-out steps.
#' @param binary_genotype_sweeps Number of genotype Gibbs sweeps per PDMI
#'   iteration.
#' @param binary_prior_gibbs_burnin Burn-in used by the binary HWE nuisance
#'   prior routines.
#' @param binary_prior_gibbs_draws Number of retained binary-prior Gibbs draws.
#' @param binary_prior_gibbs_thin Thinning interval for binary-prior Gibbs
#'   draws.
#' @param pdmi_debug If `TRUE`, retain short sampler diagnostic traces.
#' @param ncores Reserved for later parallel completed-data fits.
#' @param ... Additional named V2.2 configuration overrides.
#'
#' @return A `pdmi_control` list.
#' @export
pdmi_control <- function(agemin = 0,
                         frailtypack_maxit = 35L,
                         gh_order = 20L,
                         theta_slice_widths = NULL,
                         theta_slice_sweeps = 1L,
                         theta_slice_m = 80L,
                         continuous_slice_step = 0.7,
                         continuous_slice_m = 80L,
                         binary_genotype_sweeps = 8L,
                         binary_prior_gibbs_burnin = 30L,
                         binary_prior_gibbs_draws = 20L,
                         binary_prior_gibbs_thin = 2L,
                         pdmi_debug = FALSE,
                         ncores = 1L,
                         ...) {
  dots <- list(...)
  out <- list(
    agemin = agemin,
    frailtypack_maxit = as.integer(frailtypack_maxit),
    gh_order = as.integer(gh_order),
    theta_slice_widths = theta_slice_widths,
    theta_slice_sweeps = as.integer(theta_slice_sweeps),
    theta_slice_m = as.integer(theta_slice_m),
    continuous_slice_step = continuous_slice_step,
    continuous_slice_m = as.integer(continuous_slice_m),
    binary_genotype_sweeps = as.integer(binary_genotype_sweeps),
    binary_prior_gibbs_burnin = as.integer(binary_prior_gibbs_burnin),
    binary_prior_gibbs_draws = as.integer(binary_prior_gibbs_draws),
    binary_prior_gibbs_thin = as.integer(binary_prior_gibbs_thin),
    pdmi_debug = isTRUE(pdmi_debug),
    ncores = as.integer(ncores)
  )
  out <- c(out, dots)
  class(out) <- "pdmi_control"
  out
}

pdmi_apply_control <- function(config, control) {
  if (is.null(control)) return(config)
  if (!inherits(control, "pdmi_control")) {
    stop("`control` must be created by pdmi_control().", call. = FALSE)
  }
  for (nm in names(control)) {
    if (identical(nm, "ncores")) next
    if (!is.null(control[[nm]])) config[[nm]] <- control[[nm]]
  }
  config
}

pdmi_apply_pen_grid <- function(config, pen_grid) {
  if (is.null(pen_grid)) return(config)
  if (!is.list(pen_grid)) {
    stop("`pen_grid` must be NULL or a list with elements age, prs, gene, and k0.",
         call. = FALSE)
  }
  if (!is.null(pen_grid$age)) config$penetrance_ages <- as.numeric(pen_grid$age)
  if (!is.null(pen_grid$ages)) config$penetrance_ages <- as.numeric(pen_grid$ages)
  if (!is.null(pen_grid$prs)) config$penetrance_prs <- as.numeric(pen_grid$prs)
  if (!is.null(pen_grid$continuous)) config$penetrance_prs <- as.numeric(pen_grid$continuous)
  if (!is.null(pen_grid$gene)) config$penetrance_gene <- as.numeric(pen_grid$gene)
  if (!is.null(pen_grid$binary)) config$penetrance_gene <- as.numeric(pen_grid$binary)
  if (!is.null(pen_grid$k0)) config$penetrance_k0 <- as.numeric(pen_grid$k0)[1]
  config
}

pdmi_configure_engine <- function(M, B, numit, prior, spec, pen_grid, control) {
  config <- v22_default_config()
  config$M_imp_pdmi <- as.integer(M)
  config$M_imp_cong <- as.integer(M)
  config$mcmc_burnin <- as.integer(B)
  config$pdmi_numit <- as.integer(B) + as.integer(numit)
  config <- pdmi_apply_prior_config(config, prior, spec)
  config <- pdmi_apply_pen_grid(config, pen_grid)
  pdmi_apply_control(config, control)
}

pdmi_progress_enabled <- function(config) {
  isTRUE(config$progress)
}

pdmi_progress_stage <- function(config, label) {
  if (pdmi_progress_enabled(config) && !is.null(label) && nzchar(label)) {
    message(label)
    try(utils::flush.console(), silent = TRUE)
  }
  invisible(NULL)
}

pdmi_progress_open <- function(config, total, title = NULL, initial_label = NULL) {
  if (!pdmi_progress_enabled(config) || !is.finite(total) || total <= 0) {
    return(NULL)
  }
  if (!is.null(title) && nzchar(title)) message(title)
  env <- new.env(parent = emptyenv())
  env$total <- as.integer(total)
  env$value <- 0L
  env$closed <- FALSE
  env$message_every <- max(1L, floor(env$total / 20L))
  env$last_label_message <- NULL
  env$bar <- utils::txtProgressBar(
    min = 0,
    max = env$total,
    initial = 0,
    style = 3,
    label = initial_label %||% ""
  )
  env
}

pdmi_progress_tick <- function(progress, label = NULL, increment = 1L) {
  if (is.null(progress) || isTRUE(progress$closed)) return(invisible(progress))
  progress$value <- min(progress$total, progress$value + as.integer(increment))
  show_label <- !is.null(label) && nzchar(label) &&
    !identical(label, progress$last_label_message) &&
    (progress$value <= 1L ||
       progress$value >= progress$total ||
       grepl("; burn-in 1/", label, fixed = TRUE) ||
       grepl("; retained 1/", label, fixed = TRUE) ||
       progress$value %% progress$message_every == 0L)
  if (show_label) {
    message("  ", label)
    progress$last_label_message <- label
    try(utils::flush.console(), silent = TRUE)
  }
  utils::setTxtProgressBar(progress$bar, progress$value, label = label %||% "")
  invisible(progress)
}

pdmi_progress_close <- function(progress) {
  if (is.null(progress) || isTRUE(progress$closed)) return(invisible(NULL))
  close(progress$bar)
  progress$closed <- TRUE
  invisible(NULL)
}

pdmi_progress_iteration_label <- function(m, M, iter, total_iter, burnin) {
  burnin <- max(0L, as.integer(burnin %||% 0L))
  total_iter <- as.integer(total_iter)
  if (iter <= burnin) {
    sprintf("imputation %d/%d; burn-in %d/%d", m, M, iter, burnin)
  } else {
    sprintf("imputation %d/%d; retained %d/%d", m, M, iter - burnin,
            max(total_iter - burnin, 0L))
  }
}
