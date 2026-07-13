## ============================================================
## frailtypack complete-covariate disease analysis.
## The final disease likelihood conditions on completed covariates and
## uses frailtypack's ascertainment-corrected correlated frailty fit.
## ============================================================

v22_fit_failure <- function(reason, method = NA_character_, diagnostics = list()) {
  omega <- setNames(rep(NA_real_, 5), v22_omega_names())
  V <- matrix(NA_real_, 5, 5, dimnames = list(v22_omega_names(), v22_omega_names()))
  list(
    omega = omega,
    vcov_omega = V,
    convergence = FALSE,
    convergence_code = NA_integer_,
    failure_reason = reason,
    method = method,
    diagnostics = diagnostics,
    fit = NULL
  )
}

v22_extract_frailtypack_omega <- function(fit) {
  c(
    log.rho = log(as.numeric(fit$shape.weib[1])),
    log.lambda = log(as.numeric(fit$scale.weib[1])),
    beta_b = unname(fit$coef["mgene"]),
    beta_c = unname(fit$coef["newx"]),
    sigma_u2 = as.numeric(fit$sigma2)
  )
}

v22_match_frailtypack_coef <- function(coef_names, candidates, label) {
  hit <- match(candidates, coef_names)
  hit <- hit[!is.na(hit)]
  if (!length(hit)) {
    stop("Could not match ", label, " in fit$coef names: ",
         paste(coef_names, collapse = ", "))
  }
  hit[[1]]
}

v22_transform_frailtypack_covariance <- function(fit, coef_map = NULL) {
  raw <- as.numeric(fit$b)
  Vraw <- tryCatch(as.matrix(fit$varHtotal), error = function(e) NULL)
  if (!is.numeric(raw) || is.null(Vraw) || !is.matrix(Vraw)) {
    stop("frailtypack fit is missing fit$b or fit$varHtotal.")
  }
  if (length(raw) != nrow(Vraw) || nrow(Vraw) != ncol(Vraw)) {
    stop("fit$b and fit$varHtotal dimensions do not agree.")
  }
  if (any(!is.finite(raw)) || any(!is.finite(Vraw))) {
    stop("fit$b or fit$varHtotal contains non-finite values.")
  }

  p <- length(fit$coef)
  np <- length(raw)
  if (p < 2L || np <= p) stop("Unexpected frailtypack parameter layout.")

  ## frailtypack 2.13 optimizes positive Weibull and frailty-variance
  ## parameters on a signed square-root scale:
  ##   rho = b_rho^2, lambda = b_lambda^2, sigma_u2 = b_sigma^2.
  ## Use the signed b values so cross-covariances are transformed correctly.
  idx_shape_raw <- 1L
  idx_scale_raw <- 2L
  idx_sigma_raw <- np - p
  idx_coef <- seq.int(np - p + 1L, np)

  if (raw[idx_shape_raw] == 0 || raw[idx_scale_raw] == 0) {
    stop("Cannot transform Weibull covariance because a raw square-root parameter is zero.")
  }

  coef_names <- names(fit$coef)
  if (is.null(coef_names) || length(coef_names) != p || any(!nzchar(coef_names))) {
    stop("fit$coef must have names for covariance transformation.")
  }
  beta_b_candidates <- c(coef_map$beta_b %||% character(0),
                         "mgene", "beta_mgene", "beta_b", "majorgene")
  beta_c_candidates <- c(coef_map$beta_c %||% character(0),
                         "newx", "beta_PRS", "beta_c", "PRS")
  idx_beta_b <- idx_coef[v22_match_frailtypack_coef(coef_names, beta_b_candidates, "beta_b")]
  idx_beta_c <- idx_coef[v22_match_frailtypack_coef(coef_names, beta_c_candidates, "beta_c")]

  J <- matrix(0, nrow = length(v22_omega_names()), ncol = np)
  rownames(J) <- v22_omega_names()
  colnames(J) <- paste0("b", seq_len(np))
  colnames(J)[c(idx_shape_raw, idx_scale_raw, idx_sigma_raw,
                idx_beta_b, idx_beta_c)] <- c("sqrt.rho", "sqrt.lambda",
                                              "sqrt.sigma_u2", "beta_b", "beta_c")

  J["log.rho", idx_shape_raw] <- 2 / raw[idx_shape_raw]
  J["log.lambda", idx_scale_raw] <- 2 / raw[idx_scale_raw]
  J["beta_b", idx_beta_b] <- 1
  J["beta_c", idx_beta_c] <- 1
  J["sigma_u2", idx_sigma_raw] <- 2 * raw[idx_sigma_raw]

  V_reported <- J %*% Vraw %*% t(J)
  V_reported <- 0.5 * (V_reported + t(V_reported))
  dimnames(V_reported) <- list(v22_omega_names(), v22_omega_names())

  list(
    theta = v22_extract_frailtypack_omega(fit)[v22_omega_names()],
    vcov = V_reported[v22_omega_names(), v22_omega_names(), drop = FALSE],
    jacobian = J
  )
}

v22_extract_frailtypack_vcov <- function(fit) {
  nm <- v22_omega_names()
  transformed <- tryCatch(v22_transform_frailtypack_covariance(fit),
                          error = function(e) NULL)
  if (is.null(transformed) || any(!is.finite(transformed$vcov))) {
    return(matrix(NA_real_, 5, 5, dimnames = list(nm, nm)))
  }
  transformed$vcov[nm, nm, drop = FALSE]
}

v22_validate_analysis_dat <- function(dat, K, context = list(), kinship_cache = NULL) {
  if (!all(c("t0", "time", "status", "mgene", "newx", "famID",
             "proband", "currentage", "indID") %in% names(dat))) {
    stop("Analysis data are missing required frailtypack columns.")
  }
  finite_diag <- v22_analysis_data_diagnostics(dat, K, context = context,
                                               kinship_cache = kinship_cache)
  if (isTRUE(finite_diag$has_bad_input)) {
    bad <- paste(finite_diag$bad_columns, collapse = ", ")
    stop("Non-finite analysis data before frailtypack; bad column(s): ", bad)
  }
  if (!v22_check_popplus_support(dat)) stop("Analysis data violate pop+ support.")
  if (!v22_kinship_cache_matches(kinship_cache, dat)) {
    K <- v22_align_K(K, dat)
    blocks <- v22_family_blocks(dat)
    for (idx in blocks) {
      Ki <- K[idx, idx, drop = FALSE]
      if (!v22_is_psd(Ki, tol = 1e-6)) stop("A family K_i is not positive semidefinite.")
      tryCatch(v22_safe_chol(Ki), error = function(e) stop("A family K_i is not positive definite."))
    }
  }
  invisible(TRUE)
}

v22_fit_frailtypack <- function(dat, K, config = v22_default_config(), init_omega = NULL,
                               method = "frailtypack", context = list()) {
  v22_require_packages(include_frailtypack = TRUE)
  dat <- as.data.frame(dat)
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  kinship_cache <- config$kinship_cache %||% NULL
  K <- v22_get_cached_K(config, K, dat)
  context <- c(list(method = method), context %||% list())
  input_diag <- v22_analysis_data_diagnostics(dat, K, context = context,
                                              kinship_cache = kinship_cache)
  check <- tryCatch(v22_validate_analysis_dat(dat, K, context = context,
                                              kinship_cache = kinship_cache),
                    error = function(e) e)
  if (inherits(check, "error")) {
    return(v22_fit_failure(conditionMessage(check), method,
                           diagnostics = list(input_diagnostics = input_diag)))
  }

  sigma_lower <- as.numeric(config$frailtypack_sigma_lower %||% 1e-8)
  init_b <- NULL
  init_theta <- NULL
  init_theta_source <- "frailtypack_default"
  if (!is.null(init_omega)) {
    init_omega <- init_omega[v22_omega_names()]
    init_b <- c(beta_b = unname(init_omega["beta_b"]),
                beta_c = unname(init_omega["beta_c"]))
    init_theta <- unname(init_omega["sigma_u2"])
    init_theta_source <- "init_omega"
    if (!is.finite(init_theta) || init_theta <= sigma_lower) {
      init_theta <- NULL
      init_theta_source <- "omitted_after_boundary_start"
    }
  }

  cluster <- survival::cluster
  Surv <- survival::Surv
  fit <- tryCatch({
    null <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
    con <- file(null, open = "wt")
    sink(con)
    sink(con, type = "message")
    on.exit({
      sink(type = "message")
      sink()
      close(con)
    }, add = TRUE)
    fit_args <- list(
      formula = Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = dat,
      hazard = "Weibull",
      RandDist = "LogN",
      print.times = FALSE,
      covMatrix1 = as.matrix(K),
      recurrentAG = TRUE,
      maxit = config$frailtypack_maxit,
      proband = dat$proband,
      currentage = dat$currentage,
      agemin = config$agemin
    )
    if (!is.null(init_b)) fit_args$init.B <- init_b
    if (!is.null(init_theta)) fit_args$init.Theta <- init_theta
    suppressWarnings(do.call(frailtypack::frailtyPenal, fit_args))
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    return(v22_fit_failure(conditionMessage(fit), method,
                           diagnostics = list(input_diagnostics = input_diag,
                                              frailtypack_error = conditionMessage(fit),
                                              init_omega = init_omega,
                                              init_theta_used = init_theta,
                                              init_theta_source = init_theta_source)))
  }

  omega <- tryCatch(v22_extract_frailtypack_omega(fit), error = function(e) NULL)
  if (is.null(omega) || any(!is.finite(omega)) || omega["sigma_u2"] <= 0) {
    return(v22_fit_failure("Non-finite frailtypack estimate.", method,
                           diagnostics = list(input_diagnostics = input_diag,
                                              raw = v22_frailtypack_raw_diagnostics(fit),
                                              init_omega = init_omega,
                                              init_theta_used = init_theta,
                                              init_theta_source = init_theta_source,
                                              omega = omega)))
  }
  if (omega["sigma_u2"] <= sigma_lower) {
    return(v22_fit_failure("Boundary frailtypack frailty variance estimate.", method,
                           diagnostics = list(input_diagnostics = input_diag,
                                              raw = v22_frailtypack_raw_diagnostics(fit),
                                              init_omega = init_omega,
                                              init_theta_used = init_theta,
                                              init_theta_source = init_theta_source,
                                              omega = omega)))
  }
  V <- v22_extract_frailtypack_vcov(fit)
  if (any(!is.finite(V))) {
    return(v22_fit_failure("Non-finite frailtypack covariance on reported scale.", method,
                           diagnostics = list(input_diagnostics = input_diag,
                                              raw = v22_frailtypack_raw_diagnostics(fit),
                                              init_omega = init_omega,
                                              init_theta_used = init_theta,
                                              init_theta_source = init_theta_source,
                                              omega = omega)))
  }
  vcov_psd_before_near <- v22_is_psd(V)
  if (!vcov_psd_before_near) V <- v22_near_psd(V)
  list(
    omega = omega[v22_omega_names()],
    vcov_omega = V[v22_omega_names(), v22_omega_names(), drop = FALSE],
    convergence = TRUE,
    convergence_code = fit$istop %||% NA_integer_,
    failure_reason = NA_character_,
    method = method,
    diagnostics = list(input_diagnostics = input_diag,
                       raw = v22_frailtypack_raw_diagnostics(fit),
                       init_omega = init_omega,
                       init_theta_used = init_theta,
                       init_theta_source = init_theta_source,
                       vcov_psd_before_near = vcov_psd_before_near),
    fit = fit
  )
}

v22_valid_omega_start <- function(omega) {
  if (is.null(omega)) return(FALSE)
  omega <- tryCatch(omega[v22_omega_names()], error = function(e) NULL)
  !is.null(omega) && length(omega) == length(v22_omega_names()) &&
    all(is.finite(omega)) && is.finite(omega["sigma_u2"]) && omega["sigma_u2"] > 0
}

v22_summarize_frailtypack_attempt <- function(label, fit) {
  raw <- fit$diagnostics$raw %||% list()
  list(
    label = label,
    convergence = isTRUE(fit$convergence),
    failure_reason = fit$failure_reason %||% NA_character_,
    omega = fit$omega %||% setNames(rep(NA_real_, 5), v22_omega_names()),
    init_theta_used = fit$diagnostics$init_theta_used %||% NA_real_,
    init_theta_source = fit$diagnostics$init_theta_source %||% NA_character_,
    raw_istop = raw$istop %||% NA_integer_,
    raw_shape_weib = raw$shape_weib %||% NA_real_,
    raw_scale_weib = raw$scale_weib %||% NA_real_,
    raw_sigma2 = raw$sigma2 %||% NA_real_,
    raw_b_all_finite = raw$b_all_finite %||% NA
  )
}

v22_fit_frailtypack_multistart <- function(dat, K, config = v22_default_config(),
                                           init_candidates = list(),
                                           method = "frailtypack",
                                           context = list()) {
  if (!length(init_candidates)) {
    init_candidates <- list(frailtypack_default = NULL)
  }
  labels <- names(init_candidates)
  if (is.null(labels)) labels <- rep("", length(init_candidates))
  labels[!nzchar(labels)] <- paste0("init_", which(!nzchar(labels)))
  names(init_candidates) <- labels

  attempts <- vector("list", length(init_candidates))
  last_fit <- NULL
  for (i in seq_along(init_candidates)) {
    label <- labels[[i]]
    init_omega <- init_candidates[[i]]
    fit <- v22_fit_frailtypack(
      dat, K, config = config, init_omega = init_omega, method = method,
      context = c(context, list(final_fit_init_label = label,
                                final_fit_init_attempt = i))
    )
    attempts[[i]] <- v22_summarize_frailtypack_attempt(label, fit)
    last_fit <- fit
    if (isTRUE(fit$convergence)) {
      fit$diagnostics$selected_init_label <- label
      fit$diagnostics$selected_init_attempt <- i
      fit$diagnostics$final_fit_init_attempts <- attempts[seq_len(i)]
      return(fit)
    }
  }

  if (is.null(last_fit)) {
    last_fit <- v22_fit_failure("No frailtypack final-fit initialization attempts were made.",
                                method = method)
  }
  last_fit$failure_reason <- paste0(
    "All frailtypack final-fit initializations failed; last failure: ",
    last_fit$failure_reason %||% NA_character_
  )
  last_fit$diagnostics$selected_init_label <- NA_character_
  last_fit$diagnostics$selected_init_attempt <- NA_integer_
  last_fit$diagnostics$final_fit_init_attempts <- attempts
  last_fit
}

v22_fit_mean_completed_initial <- function(dat, K, missing_type, config) {
  tmp <- dat
  if (missing_type == "continuous" && anyNA(tmp$newx)) {
    obs_mean <- mean(tmp$newx, na.rm = TRUE)
    if (!is.finite(obs_mean)) obs_mean <- 0
    tmp$newx[is.na(tmp$newx)] <- obs_mean
  }
  if (missing_type == "binary" && anyNA(tmp$mgene)) {
    obs_p <- mean(tmp$mgene[as.numeric(tmp$proband) != 1], na.rm = TRUE)
    fill <- as.integer(is.finite(obs_p) && obs_p >= 0.5)
    tmp$mgene[is.na(tmp$mgene)] <- fill
    tmp$mgene[as.numeric(tmp$proband) == 1] <- 1L
  }
  v22_fit_frailtypack(tmp, K, config, method = "initial_mean_completed")
}

v22_pool_rubin_omega <- function(fits, M_imp, require_all = TRUE) {
  ok <- vapply(fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  if (isTRUE(require_all) && any(!ok)) {
    stop("PDMI pooling requires all completed-data frailtypack fits to converge.")
  }
  fits <- fits[ok]
  if (!length(fits)) stop("No successful completed-data fits for Rubin pooling.")
  omega_mat <- do.call(rbind, lapply(fits, function(x) x$omega[v22_omega_names()]))
  M <- nrow(omega_mat)
  qbar <- colMeans(omega_mat)
  W <- Reduce("+", lapply(fits, function(x) x$vcov_omega[v22_omega_names(), v22_omega_names(), drop = FALSE])) / M
  B <- if (M > 1L) stats::cov(omega_mat) else matrix(0, 5, 5)
  dimnames(B) <- list(v22_omega_names(), v22_omega_names())
  T <- W + (1 + 1 / M_imp) * B
  list(omega = qbar, vcov_omega = v22_near_psd(T), W = W, B = B,
       Ubar = W, T_PDMI = T, M_success = M)
}

v22_draw_omega_posterior <- function(disease_fit, config = v22_default_config()) {
  if (!isTRUE(disease_fit$convergence)) {
    stop("Cannot draw disease parameter from a failed frailtypack fit.")
  }
  omega <- disease_fit$omega[v22_omega_names()]
  V <- disease_fit$vcov_omega[v22_omega_names(), v22_omega_names(), drop = FALSE]
  if (!is.finite(omega["sigma_u2"]) || omega["sigma_u2"] <= 0) {
    stop("Cannot draw disease parameter because sigma_u2 is not positive.")
  }
  phi_names <- c("log.rho", "log.lambda", "beta_b", "beta_c", "log.sigma_u2")
  phi <- c(omega[c("log.rho", "log.lambda", "beta_b", "beta_c")],
           log.sigma_u2 = log(unname(omega["sigma_u2"])))
  J <- diag(5)
  rownames(J) <- phi_names
  colnames(J) <- v22_omega_names()
  J["log.sigma_u2", "sigma_u2"] <- 1 / unname(omega["sigma_u2"])
  Vphi <- J %*% V %*% t(J)
  Vphi <- v22_near_psd(Vphi)
  draw <- v22_rmvnorm_cov(phi, Vphi)
  names(draw) <- phi_names
  omega_draw <- c(log.rho = unname(draw["log.rho"]),
                  log.lambda = unname(draw["log.lambda"]),
                  beta_b = unname(draw["beta_b"]),
                  beta_c = unname(draw["beta_c"]),
                  sigma_u2 = exp(unname(draw["log.sigma_u2"])))[v22_omega_names()]
  attr(omega_draw, "draw_diagnostics") <- list(
    phi_mean = phi,
    phi_draw = draw,
    Vphi_diag = diag(Vphi),
    log_sigma_u2_mean = unname(phi["log.sigma_u2"]),
    log_sigma_u2_draw = unname(draw["log.sigma_u2"])
  )
  omega_draw
}

v22_penetrance_from_fit <- function(omega, vcov_omega, config = v22_default_config()) {
  grid <- expand.grid(
    age = config$penetrance_ages,
    prs = config$penetrance_prs,
    gene = config$penetrance_gene,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$estimate <- mapply(v22_penetrance, grid$age, grid$prs, grid$gene,
                          MoreArgs = list(omega = omega, k0 = config$penetrance_k0,
                                          agemin = config$agemin,
                                          gh_order = config$gh_order))
  grid$se <- mapply(function(age, prs, gene) {
    g <- v22_penetrance_gradient_omega(age, prs, gene, omega,
                                      k0 = config$penetrance_k0,
                                      agemin = config$agemin,
                                      gh_order = config$gh_order)
    sqrt(max(as.numeric(t(g) %*% vcov_omega[v22_omega_names(), v22_omega_names()] %*% g), 0))
  }, grid$age, grid$prs, grid$gene)
  grid
}

v22_penetrance_from_imputed_fits <- function(fits, pooled_vcov_omega = NULL, M_imp,
                                            config = v22_default_config()) {
  ok <- vapply(fits, function(x) isTRUE(x$convergence) && all(is.finite(x$omega)), logical(1))
  if (any(!ok)) stop("PDMI penetrance pooling requires all completed-data fits to converge.")
  fits <- fits[ok]
  if (!length(fits)) stop("No successful fits for penetrance.")
  grid0 <- expand.grid(
    age = config$penetrance_ages,
    prs = config$penetrance_prs,
    gene = config$penetrance_gene,
    KEEP.OUT.ATTRS = FALSE
  )
  rows <- lapply(seq_len(nrow(grid0)), function(r) {
    qmat <- do.call(rbind, lapply(fits, function(f) {
      v22_penetrance(grid0$age[r], grid0$prs[r], grid0$gene[r], f$omega,
                    k0 = config$penetrance_k0, agemin = config$agemin,
                    gh_order = config$gh_order)
    }))
    est <- mean(qmat)
    Bq <- if (length(qmat) > 1L) stats::var(as.numeric(qmat)) else 0
    Uq <- vapply(fits, function(f) {
      v22_penetrance_gradient_omega(grid0$age[r], grid0$prs[r], grid0$gene[r],
                                   f$omega, k0 = config$penetrance_k0,
                                   agemin = config$agemin,
                                   gh_order = config$gh_order)
    }, numeric(5))
    U_within <- vapply(seq_along(fits), function(i) {
      g <- Uq[, i]
      V <- fits[[i]]$vcov_omega[v22_omega_names(), v22_omega_names(), drop = FALSE]
      as.numeric(t(g) %*% V %*% g)
    }, numeric(1))
    var <- mean(U_within) + (1 + 1 / M_imp) * Bq
    data.frame(grid0[r, , drop = FALSE], estimate = est, se = sqrt(max(var, 0)))
  })
  do.call(rbind, rows)
}
