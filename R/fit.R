#' Fit proposed PDMI for an ascertainment-corrected frailty model
#'
#' @description
#' `pdmi_frailty()` is the main package entry point. Users provide the
#' frailtypack-style analysis formula, data, kinship matrix, missing covariates,
#' and working covariate prior. The function derives the technical-note
#' congenial imputation kernel automatically and runs the proposed V2.2
#' posterior-draw MI engine.
#'
#' @param formula Analysis formula. The validated engine requires
#'   `Surv(t0, time, status) ~ binary + continuous + cluster(family)` or
#'   `Surv(time, status) ~ binary + continuous + cluster(family)`.
#' @param data Analysis data frame.
#' @param kinship Square kinship matrix for analysis records. Row/column names
#'   should match `id`; otherwise dimensions must match `data`.
#' @param impute Character vector or list naming the covariate(s) subject to
#'   missingness. Use `list(continuous = "newx")`,
#'   `list(binary = "mgene")`, or both for explicit roles. If a named
#'   covariate has no missing values in `data`, it is treated as observed with
#'   a warning; covariates with missing values must still be listed.
#' @param prior Working covariate prior created by `pdmi_prior()`.
#' @param M Number of imputations.
#' @param B Number of sampler burn-in updates before the retained update
#'   iterations.
#' @param numit Number of retained within-imputation update iterations. The
#'   V2.2 engine returns the final state after `B + numit` updates.
#' @param report One of `"both"`, `"parameters"`, or `"penetrance"`. This
#'   controls default printing; fitted objects retain all engine output.
#' @param pen_grid Optional list overriding the penetrance grid. Recognized
#'   elements are `age`/`ages`, `prs`/`continuous`, `gene`/`binary`, and `k0`.
#' @param pen_ci If `TRUE`, `pen_summary()` includes confidence intervals.
#' @param conf.level Confidence level for parameter and penetrance intervals.
#' @param id Individual id column.
#' @param proband Proband indicator column.
#' @param currentage Current age column used by ascertainment correction.
#' @param pedigree Optional pedigree data. For binary or joint missingness it
#'   should contain ids, family ids, and preferably `fatherID`/`motherID`.
#' @param seed Optional random seed.
#' @param control Control list from `pdmi_control()`.
#' @param keep_imputations If `TRUE`, retain completed imputed data sets and
#'   sampler arrays in the fitted object.
#' @param progress If `TRUE`, print progress bars for the PDMI sampler and
#'   final completed-data frailtypack fits. Defaults to `interactive()`.
#'
#' @return A `pdmi_frailty` object with pooled parameter estimates,
#'   penetrance estimates, diagnostics, and the derived imputation model.
#' @export
#'
#' @examples
#' \dontrun{
#' data("ikun_example_joint_mar20", package = "ikunfrailty")
#' data("ikun_example_kinship", package = "ikunfrailty")
#' data("ikun_example_pedigree", package = "ikunfrailty")
#'
#' prior <- pdmi_prior(
#'   continuous = list(newx = normal_kinship(~ mgene, covariance = "kinship+iid")),
#'   binary = list(mgene = carrier_hwe(q = "estimate"))
#' )
#'
#' fit <- pdmi_frailty(
#'   survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
#'   data = ikun_example_joint_mar20,
#'   kinship = ikun_example_kinship,
#'   impute = list(continuous = "newx", binary = "mgene"),
#'   prior = prior,
#'   M = 20,
#'   B = 50,
#'   numit = 10,
#'   pedigree = ikun_example_pedigree,
#'   progress = TRUE
#' )
#' summary(fit)
#' pen_summary(fit)
#' pen_plot(fit)
#' }
pdmi_frailty <- function(formula,
                         data,
                         kinship,
                         impute,
                         prior,
                         M = 20,
                         B = 50,
                         numit = 10,
                         report = c("both", "parameters", "penetrance"),
                         pen_grid = NULL,
                         pen_ci = TRUE,
                         conf.level = 0.95,
                         id = "indID",
                         proband = "proband",
                         currentage = "currentage",
                         pedigree = NULL,
                         seed = NULL,
                         control = pdmi_control(),
                         keep_imputations = FALSE,
                         progress = interactive()) {
  call <- match.call()
  report <- match.arg(report)
  show_progress <- isTRUE(progress)
  if (!inherits(prior, "pdmi_prior")) {
    stop("`prior` must be created by pdmi_prior().", call. = FALSE)
  }
  if (!is.numeric(M) || length(M) != 1L || M < 1) {
    stop("`M` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 0) {
    stop("`B` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(numit) || length(numit) != 1L || numit < 1) {
    stop("`numit` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a single number in (0, 1).", call. = FALSE)
  }

  progress_config <- list(progress = show_progress)
  pdmi_progress_stage(progress_config, "PDMI setup: parsing formula and checking inputs")
  parsed <- pdmi_parse_model_spec(
    formula = formula,
    data = data,
    kinship = kinship,
    impute = impute,
    prior = prior,
    id = id,
    proband = proband,
    currentage = currentage,
    pedigree = pedigree
  )
  spec <- parsed$spec
  prior_version <- pdmi_resolve_prior_version(prior, spec)
  pdmi_progress_stage(
    progress_config,
    sprintf(
      "PDMI setup: missing type = %s, prior = %s, M = %d, burn-in = %d, retained = %d",
      spec$missing_type, prior_version, as.integer(M), as.integer(B), as.integer(numit)
    )
  )
  config <- pdmi_configure_engine(
    M = M,
    B = B,
    numit = numit,
    prior = prior,
    spec = spec,
    pen_grid = pen_grid,
    control = control
  )
  config$progress <- show_progress
  pdmi_progress_stage(config, "PDMI setup: deriving congenial imputation model")
  model <- pdmi_build_kernel_model(spec, prior, prior_version,
                                   B = as.integer(B),
                                   numit = as.integer(numit))

  pdmi_progress_stage(config, "PDMI setup: launching V2.2 exact-slice PDMI engine")
  result <- switch(
    spec$missing_type,
    continuous = v22_run_continuous_pdmi(parsed$dat, parsed$K,
                                         prior_version = prior_version,
                                         config = config,
                                         seed = seed),
    binary = v22_run_binary_pdmi(parsed$dat, parsed$K,
                                 prior_version = prior_version,
                                 pedigree_dat = parsed$pedigree,
                                 config = config,
                                 seed = seed),
    joint = v22_run_joint_pdmi(parsed$dat, parsed$K,
                               prior_version = prior_version,
                               pedigree_dat = parsed$pedigree,
                               config = config,
                               seed = seed)
  )

  if (!isTRUE(keep_imputations)) {
    result$imputations <- NULL
    result$impDatasets <- NULL
    result$coefIter <- NULL
  }

  out <- list(
    call = call,
    formula = stats::as.formula(formula),
    report = report,
    pen_ci = isTRUE(pen_ci),
    conf.level = conf.level,
    M = as.integer(M),
    B = as.integer(B),
    numit = as.integer(numit),
    seed = seed,
    prior = prior,
    prior_version = prior_version,
    missing_type = spec$missing_type,
    spec = spec,
    config = config,
    model = model,
    result = result,
    convergence = isTRUE(result$convergence),
    failure_reason = result$failure_reason %||% NA_character_,
    diagnostics = result$diagnostics %||% list(),
    keep_imputations = isTRUE(keep_imputations)
  )
  class(out) <- "pdmi_frailty"
  out
}

pdmi_result_omega <- function(x) {
  res <- x$result %||% x
  omega <- res$pooled$omega %||% res$omega
  if (is.null(omega)) omega <- stats::setNames(rep(NA_real_, length(v22_omega_names())),
                                               v22_omega_names())
  omega[v22_omega_names()]
}

pdmi_result_vcov <- function(x) {
  res <- x$result %||% x
  V <- res$pooled$vcov_omega %||% res$vcov_omega
  if (is.null(V)) {
    V <- matrix(NA_real_, length(v22_omega_names()), length(v22_omega_names()),
                dimnames = list(v22_omega_names(), v22_omega_names()))
  }
  V[v22_omega_names(), v22_omega_names(), drop = FALSE]
}

pdmi_parameter_table <- function(x, conf.level = x$conf.level %||% 0.95) {
  omega <- pdmi_result_omega(x)
  V <- pdmi_result_vcov(x)
  se <- sqrt(pmax(diag(V), 0))
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  data.frame(
    parameter = names(omega),
    estimate = as.numeric(omega),
    se = as.numeric(se[names(omega)]),
    lower = as.numeric(omega - z * se[names(omega)]),
    upper = as.numeric(omega + z * se[names(omega)]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
