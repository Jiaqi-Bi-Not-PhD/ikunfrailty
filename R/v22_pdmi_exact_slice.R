## ============================================================
## V2.2 exact-slice PDMI sampler overrides.
##
## The disease-parameter draw is made on
## (log.rho, log.lambda, beta_b, beta_c, sigma_u2), sampling sigma_u2
## directly on the positive half-line.  The implementation uses an
## augmented exact-invariant Gibbs/slice target with latent frailties U:
##   p(theta, U | Y, X, A=1) proportional to
##   prior(theta) prod_i L_i(Y_i | X_i, U_i; theta) phi(U_i; sigma_u2 K_i)
##   / alpha_i^F(X_i; theta).
## Marginalizing U gives the frailtypack completed-data disease posterior
## induced by the same fixed likelihood convention.  This avoids the V2.1
## normal draw from fit$b/fit$varHtotal inside the imputation sampler.
## ============================================================

v22_exact_prior_version <- function() "exact_slice_weak_proper_v1"

v22_log_prior_omega_exact <- function(omega, config = v22_default_config()) {
  omega <- omega[v22_omega_names()]
  if (any(!is.finite(omega)) || omega["sigma_u2"] <= 0) return(-Inf)
  sd_log_baseline <- as.numeric(config$theta_prior_sd_log_baseline %||% 10)
  sd_beta <- as.numeric(config$theta_prior_sd_beta %||% 10)
  tau_shape <- as.numeric(config$theta_prior_tau_shape %||% 1.1)
  tau_rate <- as.numeric(config$theta_prior_tau_rate %||% 0.2)
  sum(stats::dnorm(omega[c("log.rho", "log.lambda")], 0, sd_log_baseline, log = TRUE)) +
    sum(stats::dnorm(omega[c("beta_b", "beta_c")], 0, sd_beta, log = TRUE)) +
    stats::dgamma(unname(omega["sigma_u2"]), shape = tau_shape, rate = tau_rate, log = TRUE)
}

v22_check_completed_covariates <- function(dat, K, context = list()) {
  diag <- v22_analysis_data_diagnostics(dat, K, context = context)
  if (isTRUE(diag$has_bad_input)) {
    v22_stop_pdmi_diagnostic(
      "V2.2 exact-slice PDMI has non-finite completed covariates.",
      list(stage = context, finite_diagnostics = diag, root_cause = "input_nonfinite")
    )
  }
  invisible(TRUE)
}

v22_family_alpha_log <- function(dat, K, omega, config = v22_default_config(),
                                 prob_idx = NULL) {
  K <- v22_align_K(K, dat)
  prob_idx <- prob_idx %||% v22_proband_indices(dat)
  if (anyNA(prob_idx)) return(-Inf)
  vals <- vapply(prob_idx, function(ii) {
    alpha <- v22_alpha_popplus(dat$currentage[ii], dat$newx[ii], dat$mgene[ii],
                               K[ii, ii], omega, config$agemin, config$gh_order)
    if (!is.finite(alpha) || alpha <= 0) return(NA_real_)
    log(alpha)
  }, numeric(1))
  if (any(!is.finite(vals))) -Inf else sum(vals)
}

v22_log_augmented_disease_posterior <- function(dat, K, omega, U,
                                                config = v22_default_config(),
                                                blocks_info = NULL,
                                                prob_idx = NULL) {
  omega <- omega[v22_omega_names()]
  if (any(!is.finite(omega)) || omega["sigma_u2"] <= 0) return(-Inf)
  if (v22_omega_is_extreme(omega)) return(-Inf)
  if (length(U) != nrow(dat) || any(!is.finite(U))) return(-Inf)
  prior <- v22_log_prior_omega_exact(omega, config)
  if (!is.finite(prior)) return(-Inf)

  th <- v22_theta_from_omega(omega)
  H <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  logh <- v22_log_h0(dat$time, omega, config$agemin)
  delta <- as.numeric(dat$status)
  eta <- th$beta_c * as.numeric(dat$newx) + th$beta_b * as.numeric(dat$mgene) + U
  if (any(!is.finite(H)) || any(!is.finite(logh)) || any(!is.finite(eta))) return(-Inf)
  log_surv <- sum(delta * (logh + eta) - H * v22_safe_exp(eta))
  if (!is.finite(log_surv)) return(-Inf)

  log_alpha <- v22_family_alpha_log(dat, K, omega, config, prob_idx = prob_idx)
  if (!is.finite(log_alpha)) return(-Inf)

  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat)
  log_u <- 0
  tau <- unname(omega["sigma_u2"])
  for (info in blocks_info) {
    ui <- U[info$idx]
    quad <- as.numeric(t(ui) %*% info$K_inv %*% ui)
    log_u <- log_u - 0.5 * (info$n * log(2 * pi * tau) + info$logdetK + quad / tau)
    if (!is.finite(log_u)) return(-Inf)
  }
  prior + log_surv + log_u - log_alpha
}

v22_log_marginal_family_laplace <- function(dat, K_info, omega,
                                            config = v22_default_config(),
                                            maxit = 40L, tol = 1e-7) {
  omega <- omega[v22_omega_names()]
  if (any(!is.finite(omega)) || omega["sigma_u2"] <= 0) return(-Inf)
  th <- v22_theta_from_omega(omega)
  idx <- K_info$idx
  H <- v22_H0_diff(dat$time[idx], dat$t0[idx], omega, config$agemin)
  logh <- v22_log_h0(dat$time[idx], omega, config$agemin)
  delta <- as.numeric(dat$status[idx])
  eta0 <- th$beta_c * as.numeric(dat$newx[idx]) + th$beta_b * as.numeric(dat$mgene[idx])
  if (any(!is.finite(H)) || any(!is.finite(logh)) || any(!is.finite(eta0))) return(-Inf)
  tau <- unname(omega["sigma_u2"])
  K_inv_tau <- K_info$K_inv / tau
  log_const <- -0.5 * (K_info$n * log(2 * pi * tau) + K_info$logdetK)
  f_u <- function(u) {
    eta <- eta0 + u
    if (any(!is.finite(eta))) return(-Inf)
    r <- H * v22_safe_exp(eta)
    val <- sum(delta * (logh + eta) - r) + log_const -
      0.5 * as.numeric(t(u) %*% K_inv_tau %*% u)
    if (is.finite(val)) val else -Inf
  }
  u <- rep(0, K_info$n)
  f_cur <- f_u(u)
  if (!is.finite(f_cur)) return(-Inf)
  P <- NULL
  for (it in seq_len(maxit)) {
    eta <- eta0 + u
    r <- H * v22_safe_exp(eta)
    grad <- delta - r - as.numeric(K_inv_tau %*% u)
    P <- diag(as.numeric(r), K_info$n, K_info$n) + K_inv_tau
    step <- tryCatch(as.numeric(v22_solve_spd(P, grad)), error = function(e) rep(0, K_info$n))
    if (!all(is.finite(step))) break
    scale <- 1
    accepted <- FALSE
    for (ls in 0:12) {
      u_new <- u + scale * step
      f_new <- f_u(u_new)
      if (is.finite(f_new) && f_new >= f_cur - 1e-10) {
        u <- u_new
        f_cur <- f_new
        accepted <- TRUE
        break
      }
      scale <- scale / 2
    }
    if (max(abs(scale * step)) < tol) break
    if (!accepted) break
  }
  eta <- eta0 + u
  r <- H * v22_safe_exp(eta)
  P <- diag(as.numeric(r), K_info$n, K_info$n) + K_inv_tau
  R <- tryCatch(v22_safe_chol(P), error = function(e) NULL)
  if (is.null(R)) return(-Inf)
  logdetP <- 2 * sum(log(diag(R)))
  out <- f_cur + 0.5 * K_info$n * log(2 * pi) - 0.5 * logdetP
  if (is.finite(out)) out else -Inf
}

v22_log_marginal_disease_posterior <- function(dat, K, omega,
                                               config = v22_default_config(),
                                               blocks_info = NULL,
                                               prob_idx = NULL) {
  omega <- omega[v22_omega_names()]
  if (any(!is.finite(omega)) || omega["sigma_u2"] <= 0) return(-Inf)
  if (v22_omega_is_extreme(omega)) return(-Inf)
  prior <- v22_log_prior_omega_exact(omega, config)
  if (!is.finite(prior)) return(-Inf)
  log_alpha <- v22_family_alpha_log(dat, K, omega, config, prob_idx = prob_idx)
  if (!is.finite(log_alpha)) return(-Inf)
  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat)
  log_m <- sum(vapply(blocks_info, v22_log_marginal_family_laplace, numeric(1),
                      dat = dat, omega = omega, config = config))
  if (!is.finite(log_m)) return(-Inf)
  prior + log_m - log_alpha
}

v22_elliptical_slice <- function(x_current, mean, chol_cov, loglik_fn,
                                 max_shrink = 200L) {
  x_current <- as.numeric(x_current)
  mean <- as.numeric(mean)
  if (!length(x_current)) return(x_current)
  cur_centered <- x_current - mean
  nu <- as.numeric(t(chol_cov) %*% stats::rnorm(length(x_current)))
  logy <- loglik_fn(x_current) + log(stats::runif(1L))
  if (!is.finite(logy)) stop("Elliptical slice initialized at non-finite log likelihood.")
  theta <- stats::runif(1L, 0, 2 * pi)
  theta_min <- theta - 2 * pi
  theta_max <- theta
  for (it in seq_len(max_shrink)) {
    proposal <- mean + cur_centered * cos(theta) + nu * sin(theta)
    lp <- loglik_fn(proposal)
    if (is.finite(lp) && lp >= logy) {
      attr(proposal, "ess_shrink") <- it - 1L
      return(proposal)
    }
    if (theta < 0) theta_min <- theta else theta_max <- theta
    theta <- stats::runif(1L, theta_min, theta_max)
  }
  stop("Elliptical slice failed to find an acceptable proposal.")
}

v22_draw_frailty_ess <- function(dat, K, omega, U_current = NULL,
                                 config = v22_default_config(),
                                 blocks_info = NULL) {
  K <- v22_align_K(K, dat)
  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat)
  if (is.null(U_current)) U_current <- rep(0, nrow(dat))
  U <- as.numeric(U_current)
  th <- v22_theta_from_omega(omega)
  H <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  delta <- as.numeric(dat$status)
  base_eta <- th$beta_c * as.numeric(dat$newx) + th$beta_b * as.numeric(dat$mgene)
  shrink <- integer(length(blocks_info))
  for (b in seq_along(blocks_info)) {
    info <- blocks_info[[b]]
    idx <- info$idx
    chol_cov <- v22_safe_chol(unname(omega["sigma_u2"]) * info$K)
    ll <- function(ui) {
      eta <- base_eta[idx] + ui
      sum(delta[idx] * ui - H[idx] * v22_safe_exp(eta))
    }
    prop <- v22_elliptical_slice(U[idx], rep(0, length(idx)), chol_cov, ll,
                                max_shrink = as.integer(config$ess_max_shrink %||% 200L))
    shrink[b] <- attr(prop, "ess_shrink") %||% NA_integer_
    U[idx] <- as.numeric(prop)
  }
  list(u = U, diagnostics = list(mean_shrink = mean(shrink, na.rm = TRUE),
                                 max_shrink = max(shrink, na.rm = TRUE)))
}

v22_draw_theta_exact_slice <- function(dat, K, omega_current, U_current = NULL,
                                       config = v22_default_config(),
                                       blocks_info = NULL,
                                       prob_idx = NULL) {
  omega <- omega_current[v22_omega_names()]
  if (any(!is.finite(omega)) || omega["sigma_u2"] <= 0 || v22_omega_is_extreme(omega)) {
    omega <- v22_actual_omega(config$sigma_u2_grid[1], config)
  }
  blocks_info <- blocks_info %||% v22_precompute_K_blocks(K, dat)
  prob_idx <- prob_idx %||% v22_proband_indices(dat)
  widths <- config$theta_slice_widths %||% c(log.rho = 0.08, log.lambda = 0.08,
                                            beta_b = 0.12, beta_c = 0.12,
                                            sigma_u2 = 0.05)
  widths <- widths[v22_omega_names()]
  widths[!is.finite(widths) | widths <= 0] <- 0.1
  sweeps <- max(1L, as.integer(config$theta_slice_sweeps %||% 1L))
  max_steps <- max(10L, as.integer(config$theta_slice_m %||% 80L))
  tau_lower <- as.numeric(config$theta_tau_lower %||% 1e-6)
  tau_upper <- as.numeric(config$theta_tau_upper %||% 5)
  bounds_l <- c(log.rho = -5, log.lambda = 0, beta_b = -10, beta_c = -10, sigma_u2 = tau_lower)
  bounds_u <- c(log.rho = 5, log.lambda = 8, beta_b = 10, beta_c = 10, sigma_u2 = tau_upper)
  logpost <- function(om) {
    v22_log_marginal_disease_posterior(dat, K, om, config,
                                       blocks_info = blocks_info, prob_idx = prob_idx)
  }
  lp0 <- logpost(omega)
  if (!is.finite(lp0)) {
    v22_stop_pdmi_diagnostic(
      "V2.2 exact disease-parameter slice initialized at non-finite log posterior.",
      list(stage = list(phase = "theta_slice_init"),
           omega = omega, root_cause = "theta_slice_nonfinite")
    )
  }
  for (sweep in seq_len(sweeps)) {
    for (nm in sample(v22_omega_names())) {
      f <- function(x) {
        om <- omega
        om[nm] <- x
        logpost(om)
      }
      omega[nm] <- v22_slice_univariate(f, omega[nm], w = widths[nm], m = max_steps,
                                        lower = bounds_l[nm], upper = bounds_u[nm])
    }
  }
  out <- omega[v22_omega_names()]
  attr(out, "exact_slice_diagnostics") <- list(
    prior_version = v22_exact_prior_version(),
    logpost = logpost(omega),
    target = "frailty_marginal_selected_likelihood_laplace_logscale_prior_no_jacobian",
    theta_scale = v22_omega_names(),
    sigma_u2_sampled_directly = TRUE,
    theta_prior_jacobian_log_rho_log_lambda = FALSE,
    auxiliary_u_used_in_theta_target = FALSE,
    no_fit_b_or_varHtotal_draw = TRUE,
    theta_slice_sweeps = sweeps
  )
  out
}

v22_draw_continuous_prior_posterior_exact <- function(dat, K, prior_current = NULL,
                                                      prior_version = c("C-O", "C-R"),
                                                      config = v22_default_config()) {
  prior_version <- match.arg(prior_version)
  if (prior_version == "C-O") {
    return(list(prior_version = "C-O", eta = numeric(0), known = TRUE,
                convergence = TRUE, posterior_draw = FALSE))
  }
  K <- v22_align_K(K, dat)
  blocks <- v22_family_blocks(dat)
  x <- as.numeric(dat$newx)
  Qall <- cbind(intercept = 1, carrier = as.numeric(dat$mgene))
  if (any(!is.finite(x)) || any(!is.finite(Qall))) stop("C-R exact prior draw received non-finite completed PRS data.")
  if (is.null(prior_current) || is.null(prior_current$eta)) {
    fit0 <- v22_fit_continuous_prior(dat, K, "C-R", config)
    eta0 <- fit0$eta
  } else {
    eta0 <- prior_current$eta
  }
  ag <- log(max(as.numeric(eta0["tau_g2"] %||% config$prs_sigma2 * 0.8), 1e-6))
  ae <- log(max(as.numeric(eta0["tau_e2"] %||% config$prs_sigma2 * 0.2), 1e-6))
  gamma <- c(gamma0 = as.numeric(eta0["gamma0"] %||% 0),
             gamma_carrier = as.numeric(eta0["gamma_carrier"] %||% 0))
  gamma[!is.finite(gamma)] <- 0
  prior_sd_gamma <- as.numeric(config$cr_gamma_prior_sd %||% 10)
  log_var_prior_mean <- log(as.numeric(config$prs_sigma2 %||% 0.1) / 2)
  log_var_prior_sd <- as.numeric(config$cr_log_var_prior_sd %||% 3)
  block_info <- lapply(blocks, function(idx) {
    Ki <- K[idx, idx, drop = FALSE]
    list(idx = idx, K = Ki, Q = Qall[idx, , drop = FALSE], x = x[idx])
  })
  make_S <- function(ag, ae, info) {
    nu_g <- exp(ag)
    nu_e <- exp(ae)
    Sigma <- nu_g * info$K + nu_e * diag(length(info$idx))
    R <- v22_safe_chol(Sigma)
    list(Sigma = Sigma, Sinv = chol2inv(R), logdet = 2 * sum(log(diag(R))))
  }
  draw_gamma <- function(ag, ae) {
    Prec <- diag(1 / prior_sd_gamma^2, 2)
    rhs <- rep(0, 2)
    for (info in block_info) {
      S <- make_S(ag, ae, info)$Sinv
      Prec <- Prec + t(info$Q) %*% S %*% info$Q
      rhs <- rhs + as.numeric(t(info$Q) %*% S %*% info$x)
    }
    V <- v22_solve_spd(Prec)
    m <- as.numeric(V %*% rhs)
    stats::setNames(v22_rmvnorm_cov(m, V), c("gamma0", "gamma_carrier"))
  }
  log_xi <- function(ag, ae, gamma) {
    if (!is.finite(ag) || !is.finite(ae) || ag < -16 || ae < -16 || ag > 3 || ae > 3) return(-Inf)
    val <- stats::dnorm(ag, log_var_prior_mean, log_var_prior_sd, log = TRUE) +
      stats::dnorm(ae, log_var_prior_mean, log_var_prior_sd, log = TRUE)
    for (info in block_info) {
      ss <- make_S(ag, ae, info)
      r <- info$x - as.numeric(info$Q %*% gamma)
      val <- val - 0.5 * (ss$logdet + as.numeric(t(r) %*% ss$Sinv %*% r))
      if (!is.finite(val)) return(-Inf)
    }
    val
  }
  gamma <- draw_gamma(ag, ae)
  ag <- v22_slice_univariate(function(t) log_xi(t, ae, gamma), ag,
                             w = as.numeric(config$cr_log_var_slice_width %||% 0.5),
                             m = as.integer(config$cr_log_var_slice_m %||% 60L),
                             lower = -16, upper = 3)
  ae <- v22_slice_univariate(function(t) log_xi(ag, t, gamma), ae,
                             w = as.numeric(config$cr_log_var_slice_width %||% 0.5),
                             m = as.integer(config$cr_log_var_slice_m %||% 60L),
                             lower = -16, upper = 3)
  eta <- c(gamma0 = unname(gamma["gamma0"]),
           gamma_carrier = unname(gamma["gamma_carrier"]),
           tau_g2 = exp(ag), tau_e2 = exp(ae))
  list(prior_version = "C-R", eta = eta, known = FALSE, convergence = TRUE,
       failure_reason = NA_character_, posterior_draw = TRUE,
       diagnostics = list(log_tau_g2 = ag, log_tau_e2 = ae,
                          prior_version = v22_exact_prior_version()))
}

v22_update_continuous_missing_ess <- function(dat, K, x, omega, U, prior_fit,
                                              config = v22_default_config()) {
  K <- v22_align_K(K, dat)
  x_obs <- as.numeric(dat$newx)
  is_mis <- is.na(x_obs)
  if (!any(is_mis)) return(x_obs)
  dat_cur <- dat
  dat_cur$newx <- x
  prior_cache <- v22_build_continuous_prior_cache(dat_cur, K, prior_fit, config)
  prob_idx <- v22_proband_indices(dat)
  th <- v22_theta_from_omega(omega)
  H <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  delta <- as.numeric(dat$status)
  z <- as.numeric(dat$mgene)
  shrink <- integer(0)
  for (b in seq_along(prior_cache)) {
    info <- prior_cache[[b]]
    idx <- info$idx
    mis <- which(is_mis[idx])
    if (!length(mis)) next
    obs <- which(!is_mis[idx])
    if (length(obs)) {
      Soo <- info$Sigma[obs, obs, drop = FALSE]
      Smo <- info$Sigma[mis, obs, drop = FALSE]
      a <- as.numeric(info$mu[mis] + Smo %*% v22_solve_spd(Soo, x[idx[obs]] - info$mu[obs]))
      B <- info$Sigma[mis, mis, drop = FALSE] - Smo %*% v22_solve_spd(Soo, t(Smo))
    } else {
      a <- info$mu[mis]
      B <- info$Sigma[mis, mis, drop = FALSE]
    }
    R <- v22_safe_chol(B)
    ll <- function(xmis) {
      x_prop <- x
      x_prop[idx[mis]] <- xmis
      lp <- sum(delta[idx[mis]] * th$beta_c * xmis -
                  H[idx[mis]] * v22_safe_exp(th$beta_c * xmis +
                                               th$beta_b * z[idx[mis]] + U[idx[mis]]))
      prob_global <- prob_idx[b]
      if (prob_global %in% idx[mis]) {
        alpha <- v22_alpha_popplus(dat$currentage[prob_global], x_prop[prob_global],
                                   z[prob_global], K[prob_global, prob_global],
                                   omega, config$agemin, config$gh_order)
        lp <- lp - log(alpha)
      }
      lp
    }
    prop <- v22_elliptical_slice(x[idx[mis]], a, R, ll,
                                max_shrink = as.integer(config$ess_max_shrink %||% 200L))
    shrink <- c(shrink, attr(prop, "ess_shrink") %||% NA_integer_)
    x[idx[mis]] <- as.numeric(prop)
  }
  x[!is_mis] <- x_obs[!is_mis]
  attr(x, "ess_diagnostics") <- list(mean_shrink = mean(shrink, na.rm = TRUE),
                                     max_shrink = max(shrink, na.rm = TRUE))
  x
}

v22_draw_binary_prior_posterior_exact <- function(prior_fit, G, ped_dat,
                                                  config = v22_default_config()) {
  if (prior_fit$prior_version == "B-O") return(prior_fit)
  blocks <- v22_family_blocks(ped_dat)
  cache <- v22_build_pedigree_cache(ped_dat)
  founder_count <- max(v22_binary_founder_count(cache), 1)
  A <- v22_binary_founder_A_count(G, cache, blocks)
  q0 <- as.numeric(config$binary_beta_prior_q0 %||% config$pm)
  n0 <- as.numeric(config$binary_beta_prior_n0 %||% 50)
  aq <- max(q0 * n0, 0) + 1
  bq <- max((1 - q0) * n0, 0) + 1
  q_lower <- as.numeric(config$binary_q_lower %||% 1e-4)
  q_upper <- as.numeric(config$binary_q_upper %||% 0.25)
  n_prob_fam <- sum(vapply(blocks, function(idx) any(as.numeric(ped_dat$proband[idx]) == 1),
                           logical(1)))
  logit <- function(q) log(q / (1 - q))
  phi0 <- logit(pmin(pmax(as.numeric(prior_fit$q %||% config$pm), q_lower), q_upper))
  logpost_phi <- function(phi) {
    q <- v22_inv_logit(phi)
    if (!is.finite(q) || q <= q_lower || q >= q_upper) return(-Inf)
    p_prob_carrier <- 1 - (1 - q)^2
    if (!is.finite(p_prob_carrier) || p_prob_carrier <= 0) return(-Inf)
    (A + aq) * log(q) + (2 * founder_count - A + bq) * log1p(-q) -
      n_prob_fam * log(p_prob_carrier)
  }
  phi <- v22_slice_univariate(logpost_phi, phi0,
                              w = as.numeric(config$binary_q_logit_slice_width %||% 0.35),
                              m = as.integer(config$binary_q_logit_slice_m %||% 80L),
                              lower = logit(q_lower), upper = logit(q_upper))
  q <- v22_inv_logit(phi)
  q <- pmin(pmax(q, q_lower), q_upper)
  out <- prior_fit
  out$q <- q
  out$eta <- c(q = q)
  out$posterior_draw <- TRUE
  out$diagnostics <- c(out$diagnostics %||% list(),
                       list(founder_A_count = A, founder_count = founder_count,
                            proband_conditioned_family_count = n_prob_fam,
                            q_logit = phi,
                            q_draw_kernel = "latent_genotype_carrier_conditioned_logit_slice",
                            prior_version = v22_exact_prior_version()))
  out
}

v22_draw_continuous_pdmi <- function(dat, K, prior_version = c("C-O", "C-R"),
                                     M = 20L, numit = 10L,
                                     config = v22_default_config(), seed = NULL,
                                     init_omega = NULL, debug_context = list()) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  pdmi_progress_stage(config, "PDMI initialization: checking kinship alignment")
  K <- v22_align_K(K, dat)
  pdmi_progress_stage(config, sprintf("PDMI initialization: fitting continuous prior (%s)", prior_version))
  base_prior <- v22_fit_continuous_prior(dat, K, prior_version, config)
  is_mis <- is.na(dat$newx)
  pdmi_progress_stage(config, "PDMI initialization: initializing missing continuous covariate")
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
  }
  pdmi_progress_stage(config, "PDMI initialization: preparing kinship update blocks")
  blocks_info <- v22_get_cached_K_blocks(config, K, dat)
  prob_idx <- v22_proband_indices(dat)
  progress <- pdmi_progress_open(
    config,
    total = M * numit,
    title = sprintf("PDMI sampler: continuous missingness (%s), M = %d", prior_version, M),
    initial_label = "starting"
  )
  on.exit(pdmi_progress_close(progress), add = TRUE)
  for (m in seq_len(M)) {
    current <- dat
    current$newx <- x0
    omega_state <- init_omega %||% v22_actual_omega(config$sigma_u2_grid[1], config)
    U_state <- rep(0, nrow(dat))
    prior_state <- base_prior
    iter_draws <- vector("list", numit)
    for (iter in seq_len(numit)) {
      stage <- c(debug_context, list(missing_type = "continuous",
                                     prior_version = prior_version,
                                     phase = "exact_slice_iteration",
                                     m = m, iter = iter))
      v22_check_completed_covariates(current, K, stage)
      omega_state <- v22_draw_theta_exact_slice(current, K, omega_state, U_state,
                                                config, blocks_info, prob_idx)
      prior_state <- if (prior_version == "C-R") {
        v22_draw_continuous_prior_posterior_exact(current, K, prior_state, "C-R", config)
      } else {
        base_prior
      }
      U_draw <- v22_draw_frailty_ess(current, K, omega_state, U_state, config, blocks_info)
      U_state <- U_draw$u
      current$newx <- v22_update_continuous_missing_ess(dat, K, current$newx,
                                                        omega_state, U_state,
                                                        prior_state, config)
      coefIter[m, v22_omega_names(), iter] <- omega_state[v22_omega_names()]
      if (prior_version == "C-R") coefIter[m, names(prior_state$eta), iter] <- prior_state$eta
      iter_draws[[iter]] <- list(omega = omega_state, prior = prior_state$eta,
                                 theta_diagnostics = attr(omega_state, "exact_slice_diagnostics"),
                                 U_diagnostics = U_draw$diagnostics,
                                 x_diagnostics = attr(current$newx, "ess_diagnostics"))
      add_trace(c(stage, list(omega = omega_state,
                              theta = attr(omega_state, "exact_slice_diagnostics"))))
      pdmi_progress_tick(
        progress,
        pdmi_progress_iteration_label(m, M, iter, numit, config$mcmc_burnin)
      )
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
                          prior_version = v22_exact_prior_version(),
                          disease_draw = "frailty_marginal_exact_slice",
                          imputation_stage_frailtypack_fits = 0L,
                          theta_target = "frailty_marginal_selected_likelihood_laplace_logscale_prior_no_jacobian",
                          trace_tail = trace))
}

v22_draw_binary_pdmi <- function(dat, K, prior_version = c("B-O", "B-R"),
                                 pedigree_dat = NULL, M = 20L, numit = 10L,
                                 config = v22_default_config(), seed = NULL,
                                 init_omega = NULL, debug_context = list()) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  if (any(is.na(dat$mgene) & as.numeric(dat$proband) == 1)) {
    v22_stop_pdmi_diagnostic(
      "Binary missingness includes proband carrier status; violates pop+ support.",
      list(stage = c(debug_context, list(phase = "preflight_popplus")))
    )
  }
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  pdmi_progress_stage(config, "PDMI initialization: checking kinship alignment")
  K <- v22_align_K(K, dat)
  pdmi_progress_stage(config, "PDMI initialization: preparing pedigree")
  ped <- v22_prepare_binary_pedigree(dat, pedigree_dat)
  ped_dat <- ped$pedigree_dat
  analysis_match <- ped$analysis_match
  pdmi_progress_stage(config, sprintf("PDMI initialization: fitting binary carrier prior (%s)", prior_version))
  base_prior <- v22_fit_binary_hwe_prior(dat, prior_version, config,
                                         pedigree_dat = pedigree_dat,
                                         eb = identical(prior_version, "B-R"))
  coef_names <- c(v22_omega_names(), if (prior_version == "B-R") "q" else character(0))
  coefIter <- array(NA_real_, dim = c(M, length(coef_names), numit),
                    dimnames = list(imputation = seq_len(M), parameter = coef_names,
                                    iteration = seq_len(numit)))
  completed <- vector("list", M)
  draws <- vector("list", M)
  pdmi_progress_stage(config, "PDMI initialization: preparing kinship and pedigree update blocks")
  blocks_info <- v22_get_cached_K_blocks(config, K, dat)
  ped_blocks <- v22_family_blocks(ped_dat)
  ped_cache <- v22_build_pedigree_cache(ped_dat)
  zobs <- as.numeric(ped_dat$mgene)
  is_prob <- as.numeric(ped_dat$proband) == 1
  states_list <- lapply(seq_len(nrow(ped_dat)), function(i) {
    v22_admissible_genotypes(zobs[i], is_prob[i])
  })
  prob_idx <- v22_proband_indices(dat)
  trace <- list()
  add_trace <- function(entry) {
    if (v22_pdmi_debug_enabled(config)) {
      trace[[length(trace) + 1L]] <<- entry
      if (length(trace) > 25L) trace <<- trace[(length(trace) - 24L):length(trace)]
    }
  }
  progress <- pdmi_progress_open(
    config,
    total = M * numit,
    title = sprintf("PDMI sampler: binary missingness (%s), M = %d", prior_version, M),
    initial_label = "starting"
  )
  on.exit(pdmi_progress_close(progress), add = TRUE)
  for (m in seq_len(M)) {
    prior_state <- base_prior
    G <- v22_initialize_binary_genotypes(ped_dat, prior_state$q)
    omega_state <- init_omega %||% v22_actual_omega(config$sigma_u2_grid[1], config)
    U_state <- rep(0, nrow(dat))
    iter_draws <- vector("list", numit)
    for (iter in seq_len(numit)) {
      z_analysis <- as.integer(G[analysis_match] >= 1L)
      current <- dat
      current$mgene <- z_analysis
      current$mgene[as.numeric(current$proband) == 1] <- 1L
      stage <- c(debug_context, list(missing_type = "binary",
                                     prior_version = prior_version,
                                     phase = "exact_slice_iteration",
                                     m = m, iter = iter))
      v22_check_completed_covariates(current, K, stage)
      omega_state <- v22_draw_theta_exact_slice(current, K, omega_state, U_state,
                                                config, blocks_info, prob_idx)
      prior_state <- if (prior_version == "B-R") {
        v22_draw_binary_prior_posterior_exact(prior_state, G, ped_dat, config)
      } else {
        base_prior
      }
      U_draw <- v22_draw_frailty_ess(current, K, omega_state, U_state, config, blocks_info)
      U_state <- U_draw$u
      H_analysis <- v22_H0_diff(dat$time, dat$t0, omega_state, config$agemin)
      arrays <- v22_make_binary_update_arrays(ped_dat, dat, analysis_match, U_state,
                                              H_analysis, omega_state)
      th <- v22_theta_from_omega(omega_state)
      for (sweep in seq_len(config$binary_genotype_sweeps)) {
        G <- v22_draw_binary_genotype_sweep(ped_dat, G, prior_state$q, omega_state, arrays,
                                            config, blocks = ped_blocks, cache = ped_cache,
                                            states_list = states_list, th = th)
      }
      coefIter[m, v22_omega_names(), iter] <- omega_state[v22_omega_names()]
      if (prior_version == "B-R") coefIter[m, "q", iter] <- prior_state$q
      iter_draws[[iter]] <- list(omega = omega_state, q = prior_state$q,
                                 theta_diagnostics = attr(omega_state, "exact_slice_diagnostics"),
                                 U_diagnostics = U_draw$diagnostics)
      add_trace(c(stage, list(omega = omega_state,
                              q = prior_state$q,
                              theta = attr(omega_state, "exact_slice_diagnostics"))))
      pdmi_progress_tick(
        progress,
        pdmi_progress_iteration_label(m, M, iter, numit, config$mcmc_burnin)
      )
    }
    z_analysis <- as.integer(G[analysis_match] >= 1L)
    di <- dat
    di$mgene <- z_analysis
    di$mgene[as.numeric(di$proband) == 1] <- 1L
    completed[[m]] <- di
    draws[[m]] <- list(mgene = z_analysis, G = G, iter = iter_draws)
  }
  z_draw_mat <- do.call(cbind, lapply(draws, function(x) as.integer(x$mgene)))
  missing_analysis <- is.na(dat$mgene)
  p_mis <- if (any(missing_analysis)) rowMeans(z_draw_mat[missing_analysis, , drop = FALSE]) else numeric(0)
  carrier_count_by_m <- colSums(z_draw_mat)
  q_last_by_m <- vapply(draws, function(x) {
    tail_q <- vapply(x$iter, function(it) as.numeric(it$q %||% NA_real_), numeric(1))
    tail_q <- tail_q[is.finite(tail_q)]
    if (length(tail_q)) tail(tail_q, 1) else NA_real_
  }, numeric(1))
  list(impDatasets = completed,
       fitList = NULL,
       coefIter = coefIter,
       completed = completed,
       draws = draws,
       diagnostics = list(
         n_missing = sum(is.na(dat$mgene)),
         missing_proband_carrier = sum(is.na(dat$mgene) & as.numeric(dat$proband) == 1),
         posterior_carrier_prob_missing = p_mis,
         posterior_carrier_prob_mean = if (length(p_mis)) mean(p_mis) else NA_real_,
         posterior_carrier_prob_max = if (length(p_mis)) max(p_mis) else NA_real_,
         carrier_count_by_m = carrier_count_by_m,
         carrier_count_summary = stats::quantile(carrier_count_by_m, probs = c(0, .25, .5, .75, 1)),
         q_last_by_m = q_last_by_m,
         prior_version = v22_exact_prior_version(),
         disease_draw = "frailty_marginal_exact_slice",
         imputation_stage_frailtypack_fits = 0L,
         theta_target = "frailty_marginal_selected_likelihood_laplace_logscale_prior_no_jacobian",
         binary_q_draw = if (prior_version == "B-R") "latent_genotype_carrier_conditioned_logit_slice" else "oracle",
         trace_tail = trace
       ))
}

v22_sampler_last_omega <- function(imp, m) {
  iter <- imp$draws[[m]]$iter %||% list()
  if (!length(iter)) return(NULL)
  omega <- iter[[length(iter)]]$omega %||% NULL
  if (v22_valid_omega_start(omega)) omega[v22_omega_names()] else NULL
}

v22_final_fit_init_candidates <- function(imp, m, scenario_omega,
                                          previous_success_omega = NULL) {
  candidates <- list()
  if (v22_valid_omega_start(scenario_omega)) {
    candidates$scenario_omega <- scenario_omega[v22_omega_names()]
  }
  sampler_last <- v22_sampler_last_omega(imp, m)
  if (v22_valid_omega_start(sampler_last)) {
    candidates$sampler_last_omega <- sampler_last[v22_omega_names()]
  }
  if (v22_valid_omega_start(previous_success_omega)) {
    candidates$previous_success_omega <- previous_success_omega[v22_omega_names()]
  }
  candidates["frailtypack_default"] <- list(NULL)
  candidates
}

v22_run_continuous_pdmi <- function(dat_with_missing, K, prior_version = c("C-O", "C-R"),
                                    config = v22_default_config(), seed = NULL) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  pdmi_progress_stage(config, "PDMI setup: aligning kinship matrix")
  K <- v22_align_K(K, dat_with_missing)
  pdmi_progress_stage(config, "PDMI setup: caching per-family kinship decompositions")
  config$kinship_cache <- v22_make_kinship_cache(K, dat_with_missing)
  omega_init <- v22_actual_omega(config$sigma_u2_grid[1], config)
  M <- as.integer(config$M_imp_pdmi %||% config$M_imp_cong)
  pdmi_progress_stage(config, "PDMI setup: starting continuous missing-data sampler")
  imp <- tryCatch(v22_draw_continuous_pdmi(dat_with_missing, K, prior_version,
                                           M = M,
                                           numit = as.integer(config$pdmi_numit %||% 10L),
                                           config = config,
                                           seed = seed,
                                           init_omega = omega_init,
                                           debug_context = list(method = paste0(prior_version, "-PDMI"),
                                                                seed = seed)),
                  error = function(e) e)
  if (inherits(imp, "error")) {
    return(list(convergence = FALSE,
                failure_reason = conditionMessage(imp),
                imputations = NULL,
                diagnostics = v22_condition_diagnostics(imp)))
  }
  disease_fits <- vector("list", length(imp$completed))
  previous_success_omega <- NULL
  fit_progress <- pdmi_progress_open(
    config,
    total = length(imp$completed),
    title = "Fitting completed-data frailty models",
    initial_label = "starting"
  )
  on.exit(pdmi_progress_close(fit_progress), add = TRUE)
  for (m in seq_along(imp$completed)) {
    fit <- v22_fit_frailtypack_multistart(
      imp$completed[[m]], K = K, config = config,
      init_candidates = v22_final_fit_init_candidates(
        imp, m, scenario_omega = omega_init,
        previous_success_omega = previous_success_omega
      ),
      method = paste0(prior_version, "-PDMI"),
      context = list(missing_type = "continuous",
                     prior_version = prior_version,
                     phase = "final_completed_fit",
                     m = m)
    )
    disease_fits[[m]] <- fit
    if (isTRUE(fit$convergence)) previous_success_omega <- fit$omega
    pdmi_progress_tick(
      fit_progress,
      sprintf("completed-data fit %d/%d", m, length(imp$completed))
    )
  }
  final_status <- v22_fit_status_table(disease_fits)
  pdmi_progress_stage(config, "PDMI output: pooling parameter estimates")
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
  pdmi_progress_stage(config, "PDMI output: computing penetrance estimates")
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
                       m_success = pool$M_success,
                       imputation_stage_frailtypack_fits = 0L,
                       initializer = "scenario_sampler_last_previous_default_multistart")
  )
}

v22_run_binary_pdmi <- function(dat_with_missing, K, prior_version = c("B-O", "B-R"),
                                pedigree_dat = NULL, config = v22_default_config(),
                                seed = NULL) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  if (any(is.na(dat_with_missing$mgene) & as.numeric(dat_with_missing$proband) == 1)) {
    return(list(convergence = FALSE,
                failure_reason = "Binary missingness includes proband carrier status; violates pop+ support."))
  }
  pdmi_progress_stage(config, "PDMI setup: aligning kinship matrix")
  K <- v22_align_K(K, dat_with_missing)
  pdmi_progress_stage(config, "PDMI setup: caching per-family kinship decompositions")
  config$kinship_cache <- v22_make_kinship_cache(K, dat_with_missing)
  omega_init <- v22_actual_omega(config$sigma_u2_grid[1], config)
  M <- as.integer(config$M_imp_pdmi %||% config$M_imp_cong)
  pdmi_progress_stage(config, "PDMI setup: starting binary missing-data sampler")
  imp <- tryCatch(v22_draw_binary_pdmi(dat_with_missing, K, prior_version, pedigree_dat,
                                       M = M,
                                       numit = as.integer(config$pdmi_numit %||% 10L),
                                       config = config,
                                       seed = seed,
                                       init_omega = omega_init,
                                       debug_context = list(method = paste0(prior_version, "-PDMI"),
                                                            seed = seed)),
                  error = function(e) e)
  if (inherits(imp, "error")) {
    return(list(convergence = FALSE,
                failure_reason = conditionMessage(imp),
                imputations = NULL,
                diagnostics = v22_condition_diagnostics(imp)))
  }
  disease_fits <- vector("list", length(imp$completed))
  previous_success_omega <- NULL
  fit_progress <- pdmi_progress_open(
    config,
    total = length(imp$completed),
    title = "Fitting completed-data frailty models",
    initial_label = "starting"
  )
  on.exit(pdmi_progress_close(fit_progress), add = TRUE)
  for (m in seq_along(imp$completed)) {
    fit <- v22_fit_frailtypack_multistart(
      imp$completed[[m]], K = K, config = config,
      init_candidates = v22_final_fit_init_candidates(
        imp, m, scenario_omega = omega_init,
        previous_success_omega = previous_success_omega
      ),
      method = paste0(prior_version, "-PDMI"),
      context = list(missing_type = "binary",
                     prior_version = prior_version,
                     phase = "final_completed_fit",
                     m = m)
    )
    disease_fits[[m]] <- fit
    if (isTRUE(fit$convergence)) previous_success_omega <- fit$omega
    pdmi_progress_tick(
      fit_progress,
      sprintf("completed-data fit %d/%d", m, length(imp$completed))
    )
  }
  final_status <- v22_fit_status_table(disease_fits)
  pdmi_progress_stage(config, "PDMI output: pooling parameter estimates")
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
  pdmi_progress_stage(config, "PDMI output: computing penetrance estimates")
  pen <- tryCatch(v22_penetrance_from_imputed_fits(disease_fits, M_imp = M,
                                                   config = config),
                  error = function(e) e)
  if (inherits(pen, "error")) {
    return(list(convergence = FALSE,
                failure_reason = conditionMessage(pen),
                disease_fits = disease_fits,
                pooled = pool,
                imputations = imp,
                diagnostics = list(sampler = imp$diagnostics,
                                   final_fit_status = final_status)))
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
                       m_success = pool$M_success,
                       imputation_stage_frailtypack_fits = 0L,
                       initializer = "scenario_sampler_last_previous_default_multistart")
  )
}
