## ============================================================
## Continuous PRS frailtypack-congenial PDMI implementation.
## C-O uses known N(0, 0.1 K_i). C-R estimates
## N(Q gamma, tau_g^2 K_i + tau_e^2 I), Q=(1, mgene).
## ============================================================

v22_fit_continuous_prior <- function(dat, K, prior_version = c("C-O", "C-R"),
                                    config = v22_default_config()) {
  prior_version <- match.arg(prior_version)
  if (prior_version == "C-O") {
    return(list(
      prior_version = "C-O",
      eta = numeric(0),
      vcov = matrix(0, 0, 0),
      known = TRUE,
      convergence = TRUE,
      failure_reason = NA_character_
    ))
  }
  K <- v22_align_K(K, dat)
  blocks <- v22_family_blocks(dat)
  x <- as.numeric(dat$newx)
  Q <- cbind(intercept = 1, carrier = as.numeric(dat$mgene))
  if (any(!is.finite(Q))) {
    fallback_names <- c("gamma0", "gamma_carrier", "tau_g2", "tau_e2")
    fallback_vcov <- diag(1e6, 4)
    dimnames(fallback_vcov) <- list(fallback_names, fallback_names)
    return(list(prior_version = "C-R",
                eta = c(gamma0 = 0, gamma_carrier = 0, tau_g2 = 0.05, tau_e2 = 0.05),
                vcov = fallback_vcov,
                known = FALSE,
                convergence = FALSE,
                failure_reason = "C-R prior design matrix contains non-finite values."))
  }

  obs_block_info <- lapply(blocks, function(idx) {
    loc <- which(is.finite(x[idx]))
    if (!length(loc)) return(NULL)
    gidx <- idx[loc]
    Ki <- 0.5 * (K[gidx, gidx, drop = FALSE] + t(K[gidx, gidx, drop = FALSE]))
    list(gidx = gidx,
         K = Ki,
         Q = Q[gidx, , drop = FALSE],
         x = x[gidx])
  })
  obs_block_info <- Filter(Negate(is.null), obs_block_info)
  n_obs <- sum(vapply(obs_block_info, function(info) length(info$x), integer(1)))
  eta_names <- c("gamma0", "gamma_carrier", "tau_g2", "tau_e2")
  fallback_vcov <- function(scale = 1e6) {
    out <- diag(scale, 4)
    dimnames(out) <- list(eta_names, eta_names)
    out
  }
  if (n_obs < 5L) {
    return(list(prior_version = "C-R",
                eta = c(gamma0 = 0, gamma_carrier = 0, tau_g2 = 0.05, tau_e2 = 0.05),
                vcov = fallback_vcov(),
                known = FALSE,
                convergence = FALSE,
                failure_reason = "Too few observed PRS values to fit C-R prior."))
  }

  obs <- is.finite(x)
  lmfit <- tryCatch(stats::lm(x[obs] ~ Q[obs, 2]), error = function(e) NULL)
  gamma0 <- if (is.null(lmfit)) mean(x[obs]) else unname(coef(lmfit)[1])
  gammab <- if (is.null(lmfit) || length(coef(lmfit)) < 2L) 0 else unname(coef(lmfit)[2])
  if (!is.finite(gamma0)) gamma0 <- mean(x[obs], na.rm = TRUE)
  if (!is.finite(gamma0)) gamma0 <- 0
  if (!is.finite(gammab)) gammab <- 0
  resid <- x[obs] - as.numeric(Q[obs, , drop = FALSE] %*% c(gamma0, gammab))
  v0 <- stats::var(resid[is.finite(resid)])
  if (!is.finite(v0) || v0 <= 1e-5) v0 <- config$prs_sigma2
  init <- c(gamma0 = gamma0, gamma_carrier = gammab,
            log_tau_g2 = log(max(v0 * 0.8, 1e-5)),
            log_tau_e2 = log(max(v0 * 0.2, 1e-5)))

  nll <- function(par) {
    if (any(!is.finite(par))) return(1e100)
    gamma <- par[1:2]
    tau_g2 <- exp(par[3])
    tau_e2 <- exp(par[4])
    if (!is.finite(tau_g2) || !is.finite(tau_e2) || tau_g2 <= 0 || tau_e2 <= 0) return(1e100)
    val <- 0
    for (info in obs_block_info) {
      n_i <- length(info$x)
      Sigma <- tau_g2 * info$K + tau_e2 * diag(n_i)
      R <- tryCatch(v22_safe_chol(Sigma), error = function(e) NULL)
      if (is.null(R) || any(!is.finite(R))) return(1e100)
      mu <- as.numeric(info$Q %*% gamma)
      r <- info$x - mu
      sol <- tryCatch(backsolve(R, forwardsolve(t(R), r)), error = function(e) NULL)
      if (is.null(sol) || any(!is.finite(sol))) return(1e100)
      val <- val + sum(log(pmax(diag(R), 1e-300))) + 0.5 * sum(r * sol)
      if (!is.finite(val) || val > 1e99) return(1e100)
    }
    if (!is.finite(val)) 1e100 else val
  }

  lower <- c(-Inf, -Inf, log(1e-7), log(1e-7))
  upper <- c( Inf,  Inf, log(10),   log(10))
  opt <- tryCatch(stats::optim(init, nll, method = "L-BFGS-B", lower = lower, upper = upper,
                               hessian = TRUE,
                               control = list(maxit = 200L, factr = 1e7)),
                  error = function(e) e)
  if (inherits(opt, "error") || is.null(opt$par) || any(!is.finite(opt$par))) {
    return(list(prior_version = "C-R",
                eta = c(gamma0 = init[1], gamma_carrier = init[2],
                        tau_g2 = exp(init[3]), tau_e2 = exp(init[4])),
                vcov = fallback_vcov(),
                known = FALSE,
                convergence = FALSE,
                failure_reason = if (inherits(opt, "error")) conditionMessage(opt) else "C-R prior optimizer returned non-finite parameters."))
  }

  p <- opt$par
  eta <- c(gamma0 = p[1], gamma_carrier = p[2],
           tau_g2 = exp(p[3]), tau_e2 = exp(p[4]))
  names(eta) <- eta_names
  cov_failure <- NA_character_
  V <- tryCatch({
    Vlog <- v22_solve_spd(opt$hessian)
    if (any(!is.finite(Vlog))) stop("non-finite inverse Hessian")
    J <- diag(c(1, 1, eta["tau_g2"], eta["tau_e2"]), 4)
    Vraw <- J %*% Vlog %*% t(J)
    dimnames(Vraw) <- list(eta_names, eta_names)
    if (any(!is.finite(Vraw))) stop("non-finite transformed nuisance covariance")
    v22_near_psd(Vraw)
  }, error = function(e) {
    cov_failure <<- conditionMessage(e)
    fallback_vcov()
  })
  dimnames(V) <- list(eta_names, eta_names)

  opt_failure <- if (opt$convergence == 0) NA_character_ else paste("optim convergence", opt$convergence)
  list(
    prior_version = "C-R",
    eta = eta,
    vcov = V,
    known = FALSE,
    convergence = opt$convergence == 0,
    failure_reason = opt_failure,
    covariance_regularized = !is.na(cov_failure),
    covariance_failure_reason = cov_failure
  )
}

v22_continuous_prior_mu_sigma <- function(dat, K, prior_fit, config = v22_default_config()) {
  if (prior_fit$prior_version == "C-O") {
    mu <- rep(0, nrow(dat))
    Sigma <- config$prs_sigma2 * K
  } else {
    Q <- cbind(intercept = 1, carrier = as.numeric(dat$mgene))
    gamma <- prior_fit$eta[c("gamma0", "gamma_carrier")]
    mu <- as.numeric(Q %*% gamma)
    Sigma <- prior_fit$eta["tau_g2"] * K + prior_fit$eta["tau_e2"] * diag(nrow(K))
  }
  list(mu = mu, Sigma = 0.5 * (Sigma + t(Sigma)))
}

v22_build_continuous_prior_cache <- function(dat, K, prior_fit, config = v22_default_config()) {
  K <- v22_align_K(K, dat)
  blocks <- v22_family_blocks(dat)
  if (prior_fit$prior_version == "C-O") {
    mu_all <- rep(0, nrow(dat))
  } else {
    Q <- cbind(intercept = 1, carrier = as.numeric(dat$mgene))
    gamma <- prior_fit$eta[c("gamma0", "gamma_carrier")]
    mu_all <- as.numeric(Q %*% gamma)
  }
  lapply(blocks, function(idx) {
    Ki <- 0.5 * (K[idx, idx, drop = FALSE] + t(K[idx, idx, drop = FALSE]))
    if (prior_fit$prior_version == "C-O") {
      Sigma_i <- config$prs_sigma2 * Ki
    } else {
      Sigma_i <- prior_fit$eta["tau_g2"] * Ki + prior_fit$eta["tau_e2"] * diag(length(idx))
    }
    Sigma_i <- 0.5 * (Sigma_i + t(Sigma_i))
    R <- v22_safe_chol(Sigma_i)
    list(idx = idx,
         mu = mu_all[idx],
         Sigma = Sigma_i,
         Q = chol2inv(R),
         K_diag = diag(Ki))
  })
}

v22_initialize_continuous_prs <- function(dat, K, prior_fit, config = v22_default_config(),
                                         prior_cache = NULL) {
  x <- as.numeric(dat$newx)
  miss <- is.na(x)
  if (!any(miss)) return(x)
  prior_cache <- prior_cache %||% v22_build_continuous_prior_cache(dat, K, prior_fit, config)
  for (info in prior_cache) {
    idx <- info$idx
    mis <- which(miss[idx])
    if (!length(mis)) next
    obs <- which(!miss[idx])
    if (!length(obs)) {
      x[idx[mis]] <- info$mu[mis]
    } else {
      Soo <- info$Sigma[obs, obs, drop = FALSE]
      Smo <- info$Sigma[mis, obs, drop = FALSE]
      adj <- Smo %*% v22_solve_spd(Soo, x[idx[obs]] - info$mu[obs])
      x[idx[mis]] <- as.numeric(info$mu[mis] + adj)
    }
  }
  x
}

v22_draw_u_vg <- function(dat, K, omega, x, config = v22_default_config(),
                         blocks_info = NULL) {
  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat)
  H <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  delta <- as.numeric(dat$status)
  th <- v22_theta_from_omega(omega)
  eta_no_u <- th$beta_c * x + th$beta_b * as.numeric(dat$mgene)
  m <- numeric(nrow(dat))
  Q_save <- vector("list", length(blocks_info))
  for (it in seq_len(config$vg_maxit)) {
    w <- H * v22_safe_exp(eta_no_u + m)
    bvec <- delta - w
    m_new <- m
    for (b in seq_along(blocks_info)) {
      idx <- blocks_info[[b]]$idx
      wj <- pmax(w[idx], 1e-12)
      Q <- blocks_info[[b]]$K_inv / th$sigma_u2 + diag(wj, length(wj))
      m_new[idx] <- as.numeric(v22_solve_spd(Q, bvec[idx]))
      Q_save[[b]] <- Q
    }
    if (sqrt(sum((m_new - m)^2)) < config$vg_tol) {
      m <- m_new
      break
    }
    m <- m_new
  }
  U <- m
  for (b in seq_along(blocks_info)) {
    idx <- blocks_info[[b]]$idx
    U[idx] <- v22_rmvnorm_precision(m[idx], Q_save[[b]])
  }
  list(u = U, mean = m)
}

v22_draw_continuous_prior_posterior <- function(prior_fit, config = v22_default_config()) {
  if (prior_fit$prior_version == "C-O") return(prior_fit)
  if (!isTRUE(prior_fit$convergence)) {
    stop("Cannot draw C-R prior parameters from a failed completed-data prior fit.")
  }
  eta <- prior_fit$eta[c("gamma0", "gamma_carrier", "tau_g2", "tau_e2")]
  if (any(!is.finite(eta)) || eta["tau_g2"] <= 0 || eta["tau_e2"] <= 0) {
    stop("C-R prior estimate is not finite on the constrained scale.")
  }
  phi_names <- c("gamma0", "gamma_carrier", "log_tau_g2", "log_tau_e2")
  phi <- c(gamma0 = unname(eta["gamma0"]),
           gamma_carrier = unname(eta["gamma_carrier"]),
           log_tau_g2 = log(unname(eta["tau_g2"])),
           log_tau_e2 = log(unname(eta["tau_e2"])))
  J <- diag(4)
  rownames(J) <- phi_names
  colnames(J) <- names(eta)
  J["log_tau_g2", "tau_g2"] <- 1 / unname(eta["tau_g2"])
  J["log_tau_e2", "tau_e2"] <- 1 / unname(eta["tau_e2"])
  Vphi <- J %*% prior_fit$vcov[names(eta), names(eta), drop = FALSE] %*% t(J)
  draw <- v22_rmvnorm_cov(phi, v22_near_psd(Vphi))
  names(draw) <- phi_names
  out <- prior_fit
  out$eta <- c(gamma0 = unname(draw["gamma0"]),
               gamma_carrier = unname(draw["gamma_carrier"]),
               tau_g2 = min(max(exp(unname(draw["log_tau_g2"])), 1e-7), 10),
               tau_e2 = min(max(exp(unname(draw["log_tau_e2"])), 1e-7), 10))
  out$posterior_draw <- TRUE
  out
}

v22_update_continuous_missing_once <- function(dat, K, x, omega, prior_fit,
                                               config = v22_default_config(),
                                               prior_cache = NULL,
                                               blocks_info = NULL) {
  K <- v22_align_K(K, dat)
  x_obs <- as.numeric(dat$newx)
  is_mis <- is.na(x_obs)
  if (!any(is_mis)) return(x_obs)
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  dat_cur <- dat
  dat_cur$newx <- x
  prior_cache <- prior_cache %||% v22_build_continuous_prior_cache(dat_cur, K, prior_fit, config)
  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat_cur)
  th <- v22_theta_from_omega(omega)
  H <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  delta <- as.numeric(dat$status)
  z <- as.numeric(dat$mgene)
  prob_idx <- v22_proband_indices(dat)
  U <- v22_draw_u_vg(dat_cur, K, omega, x, config, blocks_info = blocks_info)$u
  for (b in seq_along(prior_cache)) {
    info <- prior_cache[[b]]
    idx <- info$idx
    mis_loc <- which(is_mis[idx])
    if (!length(mis_loc)) next
    Qprior <- info$Q
    mu_i <- info$mu
    prob_global <- prob_idx[b]
    for (loc in mis_loc) {
      gidx <- idx[loc]
      others <- setdiff(seq_along(idx), loc)
      qkk <- Qprior[loc, loc]
      mu_cond <- mu_i[loc]
      if (length(others)) {
        mu_cond <- mu_cond -
          as.numeric(Qprior[loc, others, drop = FALSE] %*% (x[idx[others]] - mu_i[others])) / qkk
      }
      var_cond <- 1 / max(qkk, 1e-12)
      include_alpha <- identical(gidx, prob_global)
      f <- function(xval) {
        lp <- -0.5 * (xval - mu_cond)^2 / var_cond +
          delta[gidx] * th$beta_c * xval -
          H[gidx] * v22_safe_exp(th$beta_c * xval + th$beta_b * z[gidx] + U[gidx])
        if (isTRUE(include_alpha)) {
          alpha <- v22_alpha_popplus(dat$currentage[gidx], xval, z[gidx], info$K_diag[loc],
                                     omega, config$agemin, config$gh_order)
          lp <- lp - log(alpha)
        }
        lp
      }
      x[gidx] <- v22_slice_univariate(f, x[gidx], w = config$continuous_slice_step,
                                      m = config$continuous_slice_m)
    }
  }
  x[!is_mis] <- x_obs[!is_mis]
  x
}

v22_draw_continuous_pdmi <- function(dat, K, prior_version = c("C-O", "C-R"),
                                     M = 20L, numit = 10L,
                                     config = v22_default_config(), seed = NULL,
                                     init_omega = NULL, debug_context = list()) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  K <- v22_align_K(K, dat)
  base_prior <- v22_fit_continuous_prior(dat, K, prior_version, config)
  if (prior_version == "C-R" && !isTRUE(base_prior$convergence)) {
    v22_stop_pdmi_diagnostic(
      base_prior$failure_reason %||% "Initial C-R prior fit failed.",
      list(stage = c(debug_context, list(phase = "initial_prior_fit")),
           prior_fit = base_prior)
    )
  }
  is_mis <- is.na(dat$newx)
  x0 <- v22_initialize_continuous_prs(dat, K, base_prior, config)
  coef_names <- c(v22_omega_names(),
                  if (prior_version == "C-R") c("gamma0", "gamma_carrier", "tau_g2", "tau_e2") else character(0))
  coefIter <- array(NA_real_, dim = c(M, length(coef_names), numit),
                    dimnames = list(imputation = seq_len(M), parameter = coef_names,
                                    iteration = seq_len(numit)))
  completed <- vector("list", M)
  draws <- vector("list", M)
  trace <- list()
  add_trace <- function(entry) {
    if (v22_pdmi_debug_enabled(config)) {
      trace[[length(trace) + 1L]] <<- entry
      if (length(trace) > 25L) trace <<- trace[(length(trace) - 24L):length(trace)]
    }
    invisible(NULL)
  }
  for (m in seq_len(M)) {
    current <- dat
    current$newx <- x0
    chain_init_omega <- init_omega
    iter_draws <- vector("list", numit)
    for (iter in seq_len(numit)) {
      stage <- c(debug_context, list(missing_type = "continuous",
                                     prior_version = prior_version,
                                     phase = "imputation_frailtypack",
                                     m = m, iter = iter))
      disease_fit <- v22_fit_frailtypack(current, K, config, init_omega = chain_init_omega,
                                         method = paste0(prior_version, "-PDMI imputation"),
                                         context = stage)
      add_trace(list(stage = stage,
                     convergence = disease_fit$convergence,
                     failure_reason = disease_fit$failure_reason,
                     omega = disease_fit$omega,
                     input = disease_fit$diagnostics$input_diagnostics %||% list()))
      if (!isTRUE(disease_fit$convergence)) {
        v22_stop_pdmi_diagnostic(
          disease_fit$failure_reason %||% "frailtypack imputation-stage fit failed.",
          list(stage = stage,
               trace_tail = trace,
               frailtypack = disease_fit$diagnostics,
               root_cause = v22_classify_pdmi_failure(disease_fit$failure_reason,
                                                       list(frailtypack = disease_fit$diagnostics)))
        )
      }
      omega_draw <- tryCatch(v22_draw_omega_posterior(disease_fit, config),
                             error = function(e) e)
      if (inherits(omega_draw, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(omega_draw),
          list(stage = v22_pdmi_stage_phase(stage, "omega_draw"),
               trace_tail = trace,
               frailtypack = disease_fit$diagnostics,
               current_omega = disease_fit$omega,
               root_cause = "omega_draw_extreme")
        )
      }
      if (v22_omega_is_extreme(omega_draw)) {
        v22_stop_pdmi_diagnostic(
          "Extreme or invalid disease-parameter posterior draw in PDMI imputation.",
          list(stage = v22_pdmi_stage_phase(stage, "omega_draw"),
               trace_tail = trace,
               frailtypack = disease_fit$diagnostics,
               current_omega = disease_fit$omega,
               omega_draw = omega_draw,
               omega_draw_diagnostics = attr(omega_draw, "draw_diagnostics"),
               root_cause = "omega_draw_extreme")
        )
      }
      prior_fit <- if (prior_version == "C-R") {
        v22_fit_continuous_prior(current, K, "C-R", config)
      } else {
        base_prior
      }
      prior_draw <- tryCatch(v22_draw_continuous_prior_posterior(prior_fit, config),
                             error = function(e) e)
      if (inherits(prior_draw, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(prior_draw),
          list(stage = c(stage, list(phase = "prior_draw")),
               trace_tail = trace,
               prior_fit = prior_fit)
        )
      }
      updated_x <- tryCatch(
        v22_update_continuous_missing_once(dat, K, current$newx, omega_draw, prior_draw, config),
        error = function(e) e
      )
      if (inherits(updated_x, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(updated_x),
          list(stage = c(stage, list(phase = "continuous_update")),
               trace_tail = trace,
               current_omega = disease_fit$omega,
               omega_draw = omega_draw,
               prior_draw = prior_draw$eta,
               root_cause = v22_classify_pdmi_failure(conditionMessage(updated_x),
                                                       list(omega_draw = omega_draw)))
        )
      }
      current$newx <- updated_x
      update_diag <- v22_analysis_data_diagnostics(
        current, K, context = c(stage, list(phase = "post_continuous_update"))
      )
      if (isTRUE(update_diag$has_bad_input)) {
        v22_stop_pdmi_diagnostic(
          "Non-finite analysis data after continuous PDMI update.",
          list(stage = c(stage, list(phase = "post_continuous_update")),
               trace_tail = trace,
               finite_diagnostics = update_diag,
               omega_draw = omega_draw,
               prior_draw = prior_draw$eta,
               root_cause = v22_classify_pdmi_failure(
                 "Non-finite analysis data after continuous PDMI update.",
                 list(finite_diagnostics = update_diag,
                      omega_draw = omega_draw)
               ))
        )
      }
      coefIter[m, v22_omega_names(), iter] <- omega_draw[v22_omega_names()]
      if (prior_version == "C-R") {
        coefIter[m, names(prior_draw$eta), iter] <- prior_draw$eta
      }
      chain_init_omega <- disease_fit$omega
      iter_draws[[iter]] <- list(omega = omega_draw, prior = prior_draw$eta)
    }
    completed[[m]] <- current
    draws[[m]] <- list(newx = current$newx, iter = iter_draws)
  }
  list(impDatasets = completed,
       fitList = NULL,
       coefIter = coefIter,
       completed = completed,
       draws = draws,
       diagnostics = list(n_missing = sum(is_mis),
                          missing_proband_prs = sum(is_mis & as.numeric(dat$proband) == 1),
                          trace_tail = trace))
}

v22_run_continuous_pdmi <- function(dat_with_missing, K, prior_version = c("C-O", "C-R"),
                                      config = v22_default_config(), seed = NULL) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  K <- v22_align_K(K, dat_with_missing)
  init_fit <- v22_fit_mean_completed_initial(dat_with_missing, K, "continuous", config)
  omega <- if (isTRUE(init_fit$convergence)) init_fit$omega else v22_actual_omega(config$sigma_u2_grid[1], config)
  M <- as.integer(config$M_imp_pdmi %||% config$M_imp_cong)
  imp <- tryCatch(v22_draw_continuous_pdmi(dat_with_missing, K, prior_version,
                                           M = M,
                                           numit = as.integer(config$pdmi_numit %||% 10L),
                                           config = config,
                                           seed = seed,
                                           init_omega = omega,
                                           debug_context = list(method = paste0(prior_version, "-PDMI"),
                                                                seed = seed)),
                  error = function(e) e)
  if (inherits(imp, "error")) {
    return(list(convergence = FALSE,
                failure_reason = conditionMessage(imp),
                imputations = NULL,
                diagnostics = v22_condition_diagnostics(imp)))
  }
  disease_fits <- lapply(seq_along(imp$completed), function(m) {
    v22_fit_frailtypack(imp$completed[[m]], K = K, config = config,
                        init_omega = omega, method = paste0(prior_version, "-PDMI"),
                        context = list(missing_type = "continuous",
                                       prior_version = prior_version,
                                       phase = "final_completed_fit",
                                       m = m))
  })
  final_status <- v22_fit_status_table(disease_fits)
  pool <- tryCatch(v22_pool_rubin_omega(disease_fits, M_imp = M, require_all = TRUE),
                   error = function(e) e)
  if (inherits(pool, "error")) {
    return(list(
      convergence = FALSE,
      failure_reason = conditionMessage(pool),
      disease_fits = disease_fits,
      imputations = imp,
      diagnostics = list(sampler = imp$diagnostics,
                         final_fit_status = final_status,
                         root_cause = v22_classify_pdmi_failure(conditionMessage(pool),
                                                               list(final_fit_status = final_status)))
    ))
  }
  pen <- tryCatch(
    v22_penetrance_from_imputed_fits(disease_fits, M_imp = M, config = config),
    error = function(e) e
  )
  if (inherits(pen, "error")) {
    return(list(
      convergence = FALSE,
      failure_reason = conditionMessage(pen),
      disease_fits = disease_fits,
      pooled = pool,
      imputations = imp,
      diagnostics = list(sampler = imp$diagnostics,
                         final_fit_status = final_status)
    ))
  }
  list(
    convergence = TRUE,
    failure_reason = NA_character_,
    pooled = pool,
    penetrance = pen,
    disease_fits = disease_fits,
    imputations = imp,
    impDatasets = imp$impDatasets,
    fitList = disease_fits,
    coefIter = imp$coefIter,
    diagnostics = list(sampler = imp$diagnostics,
                       final_fit_status = final_status,
                       m_success = pool$M_success)
  )
}
