#' Construct a PDMI working covariate prior
#'
#' @description
#' `pdmi_prior()` collects the user-specified working prior
#' \eqn{g_i^A(X_i \mid P_i;\eta_X)}. The disease likelihood,
#' ascertainment denominator, support factor, and frailty sampler are derived
#' from the analysis formula by `pdmi_frailty()`.
#'
#' @param continuous Named list of continuous covariate priors, typically
#'   `list(newx = normal_kinship(~ mgene, covariance = "kinship+iid"))`.
#' @param binary Named list of binary covariate priors, typically
#'   `list(mgene = carrier_hwe(q = "estimate"))`.
#'
#' @return A `pdmi_prior` object.
#' @export
pdmi_prior <- function(continuous = list(), binary = list()) {
  pdmi_check_named_prior_list(continuous, "continuous")
  pdmi_check_named_prior_list(binary, "binary")
  out <- list(continuous = continuous, binary = binary)
  class(out) <- "pdmi_prior"
  out
}

pdmi_check_named_prior_list <- function(x, label) {
  if (!is.list(x)) stop("`", label, "` must be a named list.", call. = FALSE)
  if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop("`", label, "` must be a named list with one name per covariate.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Normal kinship prior for a continuous PRS covariate
#'
#' @description
#' Describes the working prior for the continuous PRS covariate. In the V2.2
#' scalar engine, an estimated `kinship+iid` prior maps to `C-R`/`J-R`; a known
#' mean-zero `kinship` prior maps to `C-O`/`J-O`.
#'
#' @param formula Mean model for the PRS prior. The validated estimated prior
#'   uses `~ mgene`.
#' @param covariance Covariance structure. `"kinship+iid"` estimates genetic
#'   and iid variance components; `"kinship"` with `estimate = FALSE` uses a
#'   known mean-zero kinship covariance.
#' @param sigma2 Known PRS variance for the oracle `"kinship"` prior. Defaults
#'   to the V2.2 value when omitted.
#' @param estimate If `TRUE`, estimate nuisance prior parameters inside PDMI.
#'   If `FALSE`, use known nuisance parameters.
#' @param tau_g2 Optional initial or known genetic variance component.
#' @param tau_e2 Optional initial or known iid variance component.
#'
#' @return A `pdmi_continuous_prior` object.
#' @export
normal_kinship <- function(formula = ~ 1,
                           covariance = c("kinship+iid", "kinship"),
                           sigma2 = NULL,
                           estimate = NULL,
                           tau_g2 = NULL,
                           tau_e2 = NULL) {
  covariance <- match.arg(covariance)
  if (is.null(estimate)) {
    estimate <- !(identical(covariance, "kinship") && !is.null(sigma2))
  }
  out <- list(
    type = "normal_kinship",
    formula = stats::as.formula(formula),
    covariance = covariance,
    sigma2 = sigma2,
    estimate = isTRUE(estimate),
    tau_g2 = tau_g2,
    tau_e2 = tau_e2
  )
  class(out) <- c("pdmi_continuous_prior", "pdmi_prior_component")
  out
}

#' Multivariate kinship prior placeholder
#'
#' @description
#' Records a future multivariate continuous prior specification. The validated
#' V2.2 engine is scalar and will error if this prior is used in
#' `pdmi_frailty()`.
#'
#' @param formula Mean formula.
#' @param covariance Covariance description.
#' @param ... Reserved options.
#'
#' @return A `pdmi_mvn_kinship_prior` object.
#' @export
mvn_kinship <- function(formula = ~ 1, covariance = "kinship+iid", ...) {
  out <- list(
    type = "mvn_kinship",
    formula = stats::as.formula(formula),
    covariance = covariance,
    options = list(...)
  )
  class(out) <- c("pdmi_mvn_kinship_prior", "pdmi_prior_component")
  out
}

#' HWE carrier prior for a binary covariate
#'
#' @description
#' Describes the working Mendelian/HWE carrier prior for the binary carrier
#' covariate. A numeric `q` maps to `B-O`/`J-O`; `q = "estimate"` maps to
#' `B-R`/`J-R`.
#'
#' @param q Allele frequency. Use a numeric value for a known prior or
#'   `"estimate"` for posterior nuisance-parameter draws.
#' @param q0 Prior/empirical-Bayes center for `q` when estimating.
#' @param n0 Prior effective sample size for `q` when estimating.
#' @param dominant Whether the observed binary covariate is dominant carrier
#'   status. The validated engine requires `TRUE`.
#' @param founder Founder model. The validated engine requires `"hwe"`.
#'
#' @return A `pdmi_binary_prior` object.
#' @export
carrier_hwe <- function(q = "estimate",
                        q0 = 0.02,
                        n0 = 50,
                        dominant = TRUE,
                        founder = c("hwe", "free")) {
  founder <- match.arg(founder)
  if (!(is.character(q) && identical(q, "estimate")) && !is.numeric(q)) {
    stop("`q` must be numeric or the string \"estimate\".", call. = FALSE)
  }
  if (is.numeric(q) && (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1)) {
    stop("Numeric `q` must be a single value in (0, 1).", call. = FALSE)
  }
  out <- list(
    type = "carrier_hwe",
    q = q,
    q0 = q0,
    n0 = n0,
    dominant = isTRUE(dominant),
    founder = founder
  )
  class(out) <- c("pdmi_binary_prior", "pdmi_prior_component")
  out
}

#' Bernoulli GLM prior placeholder
#'
#' @description
#' Records a future binary working-prior specification. The validated V2.2
#' engine supports the Mendelian/HWE carrier prior through `carrier_hwe()` and
#' will error if a Bernoulli GLM prior is used in `pdmi_frailty()`.
#'
#' @param formula Bernoulli GLM formula.
#' @param ... Reserved options.
#'
#' @return A `pdmi_bernoulli_glm_prior` object.
#' @export
bernoulli_glm <- function(formula, ...) {
  out <- list(type = "bernoulli_glm",
              formula = stats::as.formula(formula),
              options = list(...))
  class(out) <- c("pdmi_bernoulli_glm_prior", "pdmi_prior_component")
  out
}

pdmi_prior_version_continuous <- function(prior_component) {
  if (!inherits(prior_component, "pdmi_continuous_prior")) {
    stop("The continuous prior must be created by normal_kinship().",
         call. = FALSE)
  }
  if (!identical(prior_component$type, "normal_kinship")) {
    stop("Only normal_kinship() is supported for the scalar V2.2 continuous prior.",
         call. = FALSE)
  }
  if (isTRUE(prior_component$estimate)) return("C-R")
  if (!identical(prior_component$covariance, "kinship")) {
    stop("Known continuous priors currently require covariance = \"kinship\".",
         call. = FALSE)
  }
  "C-O"
}

pdmi_prior_version_binary <- function(prior_component) {
  if (!inherits(prior_component, "pdmi_binary_prior")) {
    stop("The binary prior must be created by carrier_hwe().", call. = FALSE)
  }
  if (!identical(prior_component$type, "carrier_hwe")) {
    stop("Only carrier_hwe() is supported for the scalar V2.2 binary prior.",
         call. = FALSE)
  }
  if (!isTRUE(prior_component$dominant)) {
    stop("The validated V2.2 engine requires dominant = TRUE.", call. = FALSE)
  }
  if (!identical(prior_component$founder, "hwe")) {
    stop("The validated V2.2 engine requires founder = \"hwe\".", call. = FALSE)
  }
  if (is.numeric(prior_component$q)) "B-O" else "B-R"
}

pdmi_resolve_prior_version <- function(prior, spec) {
  miss <- spec$missing_type
  c_prior <- prior$continuous[[spec$continuous]]
  b_prior <- prior$binary[[spec$binary]]
  if (miss %in% c("continuous", "joint") && is.null(c_prior)) {
    stop("`prior$continuous` must include `", spec$continuous, "`.", call. = FALSE)
  }
  if (miss %in% c("binary", "joint") && is.null(b_prior)) {
    stop("`prior$binary` must include `", spec$binary, "`.", call. = FALSE)
  }
  if (miss == "continuous") return(pdmi_prior_version_continuous(c_prior))
  if (miss == "binary") return(pdmi_prior_version_binary(b_prior))

  cv <- pdmi_prior_version_continuous(c_prior)
  bv <- pdmi_prior_version_binary(b_prior)
  if (identical(cv, "C-O") && identical(bv, "B-O")) return("J-O")
  if (identical(cv, "C-R") && identical(bv, "B-R")) return("J-R")
  stop("Joint missingness currently supports matched oracle priors (C-O+B-O) ",
       "or matched estimated priors (C-R+B-R), not mixed prior versions.",
       call. = FALSE)
}

pdmi_apply_prior_config <- function(config, prior, spec) {
  c_prior <- prior$continuous[[spec$continuous]]
  b_prior <- prior$binary[[spec$binary]]

  if (!is.null(c_prior)) {
    if (!is.null(c_prior$sigma2)) config$prs_sigma2 <- as.numeric(c_prior$sigma2)[1]
    if (!is.null(c_prior$tau_g2)) config$prs_sigma2 <- sum(c(as.numeric(c_prior$tau_g2)[1],
                                                             as.numeric(c_prior$tau_e2 %||% 0)[1]))
  }
  if (!is.null(b_prior)) {
    if (is.numeric(b_prior$q)) config$pm <- as.numeric(b_prior$q)[1]
    config$binary_eb_q0 <- as.numeric(b_prior$q0)[1]
    config$binary_eb_n0 <- as.numeric(b_prior$n0)[1]
    config$binary_beta_prior_q0 <- as.numeric(b_prior$q0)[1]
    config$binary_beta_prior_n0 <- as.numeric(b_prior$n0)[1]
  }
  config
}
