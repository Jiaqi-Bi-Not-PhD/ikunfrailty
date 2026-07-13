## ============================================================
## Joint continuous-PRS / binary-carrier exact-slice PDMI.
## J-O-PDMI uses oracle PRS and carrier priors. J-R-PDMI uses the
## realistic PRS prior and estimated Mendelian allele-frequency prior.
## ============================================================

v22_joint_pdmi_versions <- function(prior_version = c("J-O", "J-R")) {
  prior_version <- match.arg(prior_version)
  if (prior_version == "J-O") {
    return(list(joint = "J-O", continuous = "C-O", binary = "B-O"))
  }
  list(joint = "J-R", continuous = "C-R", binary = "B-R")
}

v22_joint_continuous_prior_log_cache <- function(dat, K, prior_fit,
                                                 config = v22_default_config(),
                                                 analysis_match = NULL) {
  if (!identical(prior_fit$prior_version, "C-R")) return(list(active = FALSE))
  K <- v22_align_K(K, dat)
  blocks <- v22_family_blocks(dat)
  gamma <- prior_fit$eta[c("gamma0", "gamma_carrier")]
  out <- lapply(blocks, function(idx) {
    Ki <- K[idx, idx, drop = FALSE]
    Sigma <- prior_fit$eta["tau_g2"] * Ki + prior_fit$eta["tau_e2"] * diag(length(idx))
    R <- v22_safe_chol(Sigma)
    list(
      active = TRUE,
      famID = as.character(dat$famID[idx[1]]),
      idx = idx,
      x = as.numeric(dat$newx[idx]),
      Sinv = chol2inv(R),
      gamma = gamma,
      analysis_match = analysis_match
    )
  })
  names(out) <- vapply(out, `[[`, character(1), "famID")
  attr(out, "active") <- TRUE
  out
}

v22_joint_continuous_prior_log_family <- function(cont_info, z_analysis) {
  if (is.null(cont_info) || !length(cont_info) || isFALSE(cont_info$active %||% FALSE)) {
    return(0)
  }
  gamma <- cont_info$gamma
  mu <- unname(gamma["gamma0"]) + unname(gamma["gamma_carrier"]) * as.numeric(z_analysis)
  r <- cont_info$x - mu
  -0.5 * as.numeric(t(r) %*% cont_info$Sinv %*% r)
}

v22_joint_local_genotype_weights <- function(states, loc, Gfam, cache, q, th,
                                             arrays, fam_global_idx,
                                             cont_info = NULL,
                                             analysis_ped_locs = integer(0)) {
  logw <- rep(-Inf, length(states))
  for (s in seq_along(states)) {
    g <- states[s]
    if (is.na(cache$father[loc]) || is.na(cache$mother[loc])) {
      prior <- v22_founder_prob(g, q)
    } else {
      prior <- v22_mendel_prob(g, Gfam[cache$mother[loc]], Gfam[cache$father[loc]])
    }
    if (!(prior > 0)) next
    lp <- log(prior)
    Gprop <- Gfam
    Gprop[loc] <- g
    if (length(cache$children[[loc]])) {
      for (ch in cache$children[[loc]]) {
        if (is.na(cache$mother[ch]) || is.na(cache$father[ch])) next
        mg <- if (cache$mother[ch] == loc) g else Gprop[cache$mother[ch]]
        fg <- if (cache$father[ch] == loc) g else Gprop[cache$father[ch]]
        child_prior <- v22_mendel_prob(Gprop[ch], mg, fg)
        if (!(child_prior > 0)) {
          lp <- -Inf
          break
        }
        lp <- lp + log(child_prior)
      }
    }
    if (!is.finite(lp)) next
    z <- as.integer(g >= 1L)
    global <- fam_global_idx[loc]
    disease <- arrays$status[global] * th$beta_b * z -
      arrays$H[global] * v22_safe_exp(arrays$base_eta[global] + th$beta_b * z)
    lp <- lp + disease
    if (length(analysis_ped_locs) &&
        !is.null(cont_info) &&
        isTRUE(cont_info$active %||% FALSE)) {
      z_analysis <- as.integer(Gprop[analysis_ped_locs] >= 1L)
      lp <- lp + v22_joint_continuous_prior_log_family(cont_info, z_analysis)
    }
    logw[s] <- lp
  }
  m <- max(logw)
  if (!is.finite(m)) return(rep(1 / length(states), length(states)))
  w <- exp(logw - m)
  if (!any(is.finite(w)) || sum(w) <= 0) return(rep(1 / length(states), length(states)))
  w / sum(w)
}

v22_draw_joint_binary_genotype_sweep <- function(ped_dat, G, q, omega, arrays,
                                                 dat_current, K, continuous_prior,
                                                 analysis_match,
                                                 config = v22_default_config(),
                                                 blocks = NULL, cache = NULL,
                                                 states_list = NULL, th = NULL) {
  blocks <- blocks %||% v22_family_blocks(ped_dat)
  cache <- cache %||% v22_build_pedigree_cache(ped_dat)
  th <- th %||% v22_theta_from_omega(omega)
  Gnew <- G
  if (is.null(states_list)) {
    zobs <- as.numeric(ped_dat$mgene)
    is_prob <- as.numeric(ped_dat$proband) == 1
    states_list <- lapply(seq_len(nrow(ped_dat)), function(i) {
      v22_admissible_genotypes(zobs[i], is_prob[i])
    })
  }
  cont_cache <- v22_joint_continuous_prior_log_cache(
    dat_current, K, continuous_prior, config, analysis_match = analysis_match
  )
  analysis_blocks <- v22_family_blocks(dat_current)
  fam_to_analysis_idx <- lapply(analysis_blocks, identity)
  names(fam_to_analysis_idx) <- names(analysis_blocks)

  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]
    fam <- as.character(ped_dat$famID[idx[1]])
    Gfam <- Gnew[idx]
    cont_info <- cont_cache[[fam]] %||% NULL
    analysis_idx <- fam_to_analysis_idx[[fam]] %||% integer(0)
    analysis_ped_locs <- if (length(analysis_idx)) {
      match(analysis_match[analysis_idx], idx)
    } else integer(0)
    analysis_ped_locs <- analysis_ped_locs[!is.na(analysis_ped_locs)]
    for (loc in seq_along(idx)) {
      states <- states_list[[idx[loc]]]
      if (length(states) == 1L) {
        Gfam[loc] <- states[1]
      } else {
        w <- v22_joint_local_genotype_weights(
          states, loc, Gfam, cache[[b]], q, th, arrays, idx,
          cont_info = cont_info, analysis_ped_locs = analysis_ped_locs
        )
        Gfam[loc] <- sample(states, 1L, prob = w)
      }
    }
    Gnew[idx] <- Gfam
  }
  Gnew
}

v22_update_joint_continuous_missing_ess <- function(dat_obs, dat_current, K, x,
                                                    omega, U, prior_fit,
                                                    config = v22_default_config()) {
  K <- v22_align_K(K, dat_current)
  x_obs <- as.numeric(dat_obs$newx)
  is_mis <- is.na(x_obs)
  if (!any(is_mis)) return(x_obs)
  dat_cur <- dat_current
  dat_cur$newx <- x
  prior_cache <- v22_build_continuous_prior_cache(dat_cur, K, prior_fit, config)
  prob_idx <- v22_proband_indices(dat_current)
  th <- v22_theta_from_omega(omega)
  H <- v22_H0_diff(dat_current$time, dat_current$t0, omega, config$agemin)
  delta <- as.numeric(dat_current$status)
  z <- as.numeric(dat_current$mgene)
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
        alpha <- v22_alpha_popplus(dat_current$currentage[prob_global],
                                   x_prop[prob_global],
                                   z[prob_global],
                                   K[prob_global, prob_global],
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

v22_initialize_joint_pdmi_state <- function(dat, K, pedigree_dat, prior_version,
                                            config = v22_default_config()) {
  versions <- v22_joint_pdmi_versions(prior_version)
  pdmi_progress_stage(config, "PDMI initialization: preparing pedigree for joint imputation")
  ped <- v22_prepare_binary_pedigree(dat, pedigree_dat)
  ped_dat <- ped$pedigree_dat
  analysis_match <- ped$analysis_match
  pdmi_progress_stage(config, sprintf("PDMI initialization: fitting binary carrier prior (%s)", versions$binary))
  binary_prior <- v22_fit_binary_hwe_prior(
    dat, versions$binary, config, pedigree_dat = pedigree_dat,
    eb = identical(versions$binary, "B-R")
  )
  pdmi_progress_stage(config, "PDMI initialization: initializing binary genotypes")
  G <- v22_initialize_binary_genotypes(ped_dat, binary_prior$q)
  dat_z <- dat
  dat_z$mgene <- as.integer(G[analysis_match] >= 1L)
  dat_z$mgene[as.numeric(dat_z$proband) == 1] <- 1L
  pdmi_progress_stage(config, sprintf("PDMI initialization: fitting continuous prior (%s)", versions$continuous))
  continuous_prior <- v22_fit_continuous_prior(dat_z, K, versions$continuous, config)
  pdmi_progress_stage(config, "PDMI initialization: initializing missing continuous covariate")
  x <- v22_initialize_continuous_prs(dat_z, K, continuous_prior, config)
  dat_z$newx <- x
  list(ped = ped,
       binary_prior = binary_prior,
       continuous_prior = continuous_prior,
       G = G,
       x = x,
       completed_dat = dat_z)
}

v22_draw_joint_pdmi <- function(dat, K, prior_version = c("J-O", "J-R"),
                                pedigree_dat = NULL, M = 20L, numit = 10L,
                                config = v22_default_config(), seed = NULL,
                                init_omega = NULL, debug_context = list()) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  if (any(is.na(dat$mgene) & as.numeric(dat$proband) == 1)) {
    v22_stop_pdmi_diagnostic(
      "Joint missingness includes proband carrier status; violates pop+ support.",
      list(stage = c(debug_context, list(phase = "preflight_popplus")))
    )
  }
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  pdmi_progress_stage(config, "PDMI initialization: checking kinship alignment")
  K <- v22_align_K(K, dat)
  versions <- v22_joint_pdmi_versions(prior_version)
  pdmi_progress_stage(config, "PDMI initialization: building joint sampler state")
  init <- v22_initialize_joint_pdmi_state(dat, K, pedigree_dat, prior_version, config)
  ped_dat <- init$ped$pedigree_dat
  analysis_match <- init$ped$analysis_match
  base_binary_prior <- init$binary_prior
  base_continuous_prior <- init$continuous_prior
  x0 <- init$x

  coef_names <- c(
    v22_omega_names(),
    if (prior_version == "J-R") c("gamma0", "gamma_carrier", "tau_g2", "tau_e2", "q") else character(0)
  )
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
    title = sprintf("PDMI sampler: joint missingness (%s), M = %d", prior_version, M),
    initial_label = "starting"
  )
  on.exit(pdmi_progress_close(progress), add = TRUE)
  for (m in seq_len(M)) {
    G <- v22_initialize_binary_genotypes(ped_dat, base_binary_prior$q)
    current <- dat
    current$mgene <- as.integer(G[analysis_match] >= 1L)
    current$mgene[as.numeric(current$proband) == 1] <- 1L
    current$newx <- x0
    omega_state <- init_omega %||% v22_actual_omega(config$sigma_u2_grid[1], config)
    U_state <- rep(0, nrow(dat))
    prior_c_state <- base_continuous_prior
    prior_b_state <- base_binary_prior
    iter_draws <- vector("list", numit)
    for (iter in seq_len(numit)) {
      current$mgene <- as.integer(G[analysis_match] >= 1L)
      current$mgene[as.numeric(current$proband) == 1] <- 1L
      stage <- c(debug_context, list(missing_type = "joint",
                                     prior_version = prior_version,
                                     phase = "joint_exact_slice_iteration",
                                     m = m, iter = iter))
      v22_check_completed_covariates(current, K, stage)
      omega_state <- v22_draw_theta_exact_slice(current, K, omega_state, U_state,
                                                config, blocks_info, prob_idx)
      prior_c_state <- if (versions$continuous == "C-R") {
        v22_draw_continuous_prior_posterior_exact(current, K, prior_c_state, "C-R", config)
      } else {
        base_continuous_prior
      }
      prior_b_state <- if (versions$binary == "B-R") {
        v22_draw_binary_prior_posterior_exact(prior_b_state, G, ped_dat, config)
      } else {
        base_binary_prior
      }
      U_draw <- v22_draw_frailty_ess(current, K, omega_state, U_state, config, blocks_info)
      U_state <- U_draw$u
      current$newx <- v22_update_joint_continuous_missing_ess(
        dat, current, K, current$newx, omega_state, U_state, prior_c_state, config
      )
      H_analysis <- v22_H0_diff(current$time, current$t0, omega_state, config$agemin)
      arrays <- v22_make_binary_update_arrays(ped_dat, current, analysis_match,
                                              U_state, H_analysis, omega_state)
      th <- v22_theta_from_omega(omega_state)
      for (sweep in seq_len(config$binary_genotype_sweeps)) {
        G <- v22_draw_joint_binary_genotype_sweep(
          ped_dat, G, prior_b_state$q, omega_state, arrays,
          dat_current = current, K = K, continuous_prior = prior_c_state,
          analysis_match = analysis_match, config = config,
          blocks = ped_blocks, cache = ped_cache, states_list = states_list, th = th
        )
      }
      current$mgene <- as.integer(G[analysis_match] >= 1L)
      current$mgene[as.numeric(current$proband) == 1] <- 1L
      coefIter[m, v22_omega_names(), iter] <- omega_state[v22_omega_names()]
      if (prior_version == "J-R") {
        coefIter[m, names(prior_c_state$eta), iter] <- prior_c_state$eta
        coefIter[m, "q", iter] <- prior_b_state$q
      }
      iter_draws[[iter]] <- list(
        omega = omega_state,
        continuous_prior = prior_c_state$eta,
        q = prior_b_state$q,
        theta_diagnostics = attr(omega_state, "exact_slice_diagnostics"),
        U_diagnostics = U_draw$diagnostics,
        x_diagnostics = attr(current$newx, "ess_diagnostics")
      )
      add_trace(c(stage, list(omega = omega_state,
                              q = prior_b_state$q,
                              theta = attr(omega_state, "exact_slice_diagnostics"))))
      pdmi_progress_tick(
        progress,
        pdmi_progress_iteration_label(m, M, iter, numit, config$mcmc_burnin)
      )
    }
    current$mgene <- as.integer(G[analysis_match] >= 1L)
    current$mgene[as.numeric(current$proband) == 1] <- 1L
    completed[[m]] <- current
    draws[[m]] <- list(newx = current$newx,
                       mgene = as.integer(current$mgene),
                       G = G,
                       iter = iter_draws)
  }
  z_draw_mat <- do.call(cbind, lapply(draws, function(x) as.integer(x$mgene)))
  x_draw_mat <- do.call(cbind, lapply(draws, function(x) as.numeric(x$newx)))
  missing_z <- is.na(dat$mgene)
  missing_x <- is.na(dat$newx)
  p_mis <- if (any(missing_z)) rowMeans(z_draw_mat[missing_z, , drop = FALSE]) else numeric(0)
  list(
    impDatasets = completed,
    fitList = NULL,
    coefIter = coefIter,
    completed = completed,
    draws = draws,
    diagnostics = list(
      n_missing_x = sum(missing_x),
      n_missing_z = sum(missing_z),
      missing_proband_prs = sum(missing_x & as.numeric(dat$proband) == 1),
      missing_proband_carrier = sum(missing_z & as.numeric(dat$proband) == 1),
      posterior_carrier_prob_missing = p_mis,
      posterior_carrier_prob_mean = if (length(p_mis)) mean(p_mis) else NA_real_,
      posterior_carrier_prob_max = if (length(p_mis)) max(p_mis) else NA_real_,
      completed_prs_missing_summary = if (any(missing_x)) v22_numeric_summary(as.numeric(x_draw_mat[missing_x, ])) else numeric(0),
      prior_version = v22_exact_prior_version(),
      disease_draw = "frailty_marginal_exact_slice",
      imputation_stage_frailtypack_fits = 0L,
      theta_target = "frailty_marginal_selected_likelihood_laplace_logscale_prior_no_jacobian",
      joint_binary_update = if (prior_version == "J-R") "includes_continuous_prior_log_density" else "oracle_independent_prs_prior",
      binary_q_draw = if (prior_version == "J-R") "latent_genotype_carrier_conditioned_logit_slice" else "oracle",
      trace_tail = trace
    )
  )
}

v22_run_joint_pdmi <- function(dat_with_missing, K, prior_version = c("J-O", "J-R"),
                               pedigree_dat = NULL, config = v22_default_config(),
                               seed = NULL) {
  prior_version <- match.arg(prior_version)
  if (!is.null(seed)) set.seed(seed)
  pdmi_progress_stage(config, "PDMI setup: aligning kinship matrix")
  K <- v22_align_K(K, dat_with_missing)
  pdmi_progress_stage(config, "PDMI setup: caching per-family kinship decompositions")
  config$kinship_cache <- v22_make_kinship_cache(K, dat_with_missing)
  omega_init <- v22_actual_omega(config$sigma_u2_grid[1], config)
  M <- as.integer(config$M_imp_pdmi %||% config$M_imp_cong)
  pdmi_progress_stage(config, "PDMI setup: starting joint missing-data sampler")
  imp <- tryCatch(v22_draw_joint_pdmi(
    dat_with_missing, K, prior_version, pedigree_dat,
    M = M,
    numit = as.integer(config$pdmi_numit %||% 10L),
    config = config,
    seed = seed,
    init_omega = omega_init,
    debug_context = list(method = paste0(prior_version, "-PDMI"), seed = seed)
  ), error = function(e) e)
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
      context = list(missing_type = "joint",
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
    return(list(convergence = FALSE,
                failure_reason = conditionMessage(pool),
                disease_fits = disease_fits,
                imputations = imp,
                diagnostics = list(sampler = imp$diagnostics,
                                   final_fit_status = final_status,
                                   root_cause = v22_classify_pdmi_failure(
                                     conditionMessage(pool),
                                     list(final_fit_status = final_status)))))
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
