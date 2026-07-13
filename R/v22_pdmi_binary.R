## ============================================================
## Binary carrier frailtypack-congenial PDMI implementation.
## B-O uses the carrier-conditioned Mendelian prior with q = 0.02 as known.
## B-R estimates the HWE allele frequency q under the pop+ proband-carrier
## support and posterior-draws q during imputation.
## ============================================================

v22_founder_prob <- function(g, q) stats::dbinom(g, size = 2L, prob = q)

v22_mendel_prob <- function(g, mother_g, father_g) {
  pmom <- mother_g / 2
  pdad <- father_g / 2
  if (g == 0L) return((1 - pmom) * (1 - pdad))
  if (g == 1L) return(pmom * (1 - pdad) + (1 - pmom) * pdad)
  if (g == 2L) return(pmom * pdad)
  0
}

v22_build_pedigree_cache <- function(ped_dat, id_col = "indID",
                                    father_col = "fatherID", mother_col = "motherID") {
  blocks <- v22_family_blocks(ped_dat)
  lapply(blocks, function(idx) {
    ids <- as.character(ped_dat[[id_col]][idx])
    father <- as.character(ped_dat[[father_col]][idx])
    mother <- as.character(ped_dat[[mother_col]][idx])
    father[father %in% c("", "0", "NA")] <- NA_character_
    mother[mother %in% c("", "0", "NA")] <- NA_character_
    fi <- match(father, ids)
    mi <- match(mother, ids)
    children <- lapply(seq_along(idx), function(k) which(fi == k | mi == k))
    list(idx = idx, ids = ids, father = fi, mother = mi, children = children)
  })
}

v22_admissible_genotypes <- function(carrier_value, is_proband = FALSE) {
  if (isTRUE(is_proband)) return(c(1L, 2L))
  if (!is.finite(carrier_value)) return(0:2)
  if (as.integer(carrier_value) == 0L) return(0L)
  c(1L, 2L)
}

v22_initialize_binary_genotypes <- function(ped_dat, q, proband_col = "proband") {
  zobs <- as.numeric(ped_dat$mgene)
  is_prob <- as.numeric(ped_dat[[proband_col]]) == 1
  zobs[is_prob] <- 1
  G <- integer(nrow(ped_dat))
  for (i in seq_len(nrow(ped_dat))) {
    states <- v22_admissible_genotypes(zobs[i], is_prob[i])
    pr <- vapply(states, v22_founder_prob, numeric(1), q = q)
    if (!any(pr > 0)) pr <- rep(1, length(states))
    G[i] <- sample(states, 1L, prob = pr)
  }
  G
}

v22_binary_prior_arrays <- function(n) {
  list(
    status = rep(0, n),
    H = rep(0, n),
    base_eta = rep(0, n)
  )
}

v22_binary_founder_count <- function(cache) {
  sum(vapply(cache, function(fam) {
    sum(is.na(fam$father) | is.na(fam$mother))
  }, numeric(1)))
}

v22_binary_founder_A_count <- function(G, cache, blocks) {
  sum(vapply(seq_along(blocks), function(b) {
    loc <- which(is.na(cache[[b]]$father) | is.na(cache[[b]]$mother))
    if (!length(loc)) return(0)
    sum(G[blocks[[b]][loc]])
  }, numeric(1)))
}

v22_sample_binary_prior_genotypes <- function(ped_dat, q, config = v22_default_config()) {
  blocks <- v22_family_blocks(ped_dat)
  cache <- v22_build_pedigree_cache(ped_dat)
  zobs <- as.numeric(ped_dat$mgene)
  is_prob <- as.numeric(ped_dat$proband) == 1
  zobs[is_prob] <- 1
  ped_work <- ped_dat
  ped_work$mgene <- zobs
  states_list <- lapply(seq_len(nrow(ped_work)), function(i) {
    v22_admissible_genotypes(zobs[i], is_prob[i])
  })
  G <- v22_initialize_binary_genotypes(ped_work, q)
  arrays <- v22_binary_prior_arrays(nrow(ped_work))
  th <- list(beta_b = 0)
  burnin <- max(0L, as.integer(config$binary_prior_gibbs_burnin %||% 30L))
  draws <- max(1L, as.integer(config$binary_prior_gibbs_draws %||% 20L))
  thin <- max(1L, as.integer(config$binary_prior_gibbs_thin %||% 2L))
  total <- burnin + draws * thin
  A <- numeric(draws)
  save_i <- 0L
  for (iter in seq_len(total)) {
    G <- v22_draw_binary_genotype_sweep(
      ped_work, G, q, omega = NULL, arrays = arrays, config = config,
      blocks = blocks, cache = cache, states_list = states_list, th = th
    )
    if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
      save_i <- save_i + 1L
      A[save_i] <- v22_binary_founder_A_count(G, cache, blocks)
    }
  }
  founder_count <- v22_binary_founder_count(cache)
  score <- A / q - (2 * founder_count - A) / (1 - q)
  hessian <- -A / q^2 - (2 * founder_count - A) / (1 - q)^2
  score_var <- if (length(score) > 1L) stats::var(score) else 0
  list(
    expected_A = mean(A),
    sd_A = stats::sd(A),
    A_draws = A,
    founder_count = founder_count,
    score_q_mean = mean(score),
    info_q = max(-mean(hessian) - score_var, 0),
    q = q
  )
}

v22_binary_proband_only_pedigree <- function(ped_dat) {
  out <- ped_dat
  out$mgene <- NA_real_
  out$mgene[as.numeric(out$proband) == 1] <- 1L
  out
}

v22_estimate_binary_q_carrier_conditioned <- function(dat, pedigree_dat = NULL,
                                                     config = v22_default_config(),
                                                     eb = FALSE) {
  ## Updated Section 6 / binary Section 3.6: estimate q under the
  ## carrier-conditioned Mendelian prior, not from ordinary HWE prevalence
  ## in the observed nonproband records.
  ped <- v22_prepare_binary_pedigree(dat, pedigree_dat)$pedigree_dat
  ped$mgene[as.numeric(ped$proband) == 1] <- 1L
  eligible <- as.numeric(dat$proband) != 1 & is.finite(dat$mgene)
  n_obs <- sum(eligible)
  if (n_obs < 5L) {
    return(list(q = config$pm, vcov = 1e6, convergence = FALSE,
                failure_reason = "Too few observed nonproband carrier statuses to estimate q.",
                diagnostics = list(n_observed_nonproband = n_obs)))
  }

  p_carrier <- mean(as.numeric(dat$mgene[eligible]))
  p_carrier <- pmin(pmax(p_carrier, 1e-8), 1 - 1e-8)
  q_naive <- 1 - sqrt(1 - p_carrier)
  q_lower <- config$binary_q_lower %||% 1e-4
  q_upper <- config$binary_q_upper %||% 0.25
  eb_q0 <- pmin(pmax(config$binary_eb_q0 %||% config$pm, q_lower), q_upper)
  eb_n0 <- if (isTRUE(eb)) max(as.numeric(config$binary_eb_n0 %||% 0), 0) else 0
  q_grid <- unique(sort(c(
    exp(seq(log(q_lower), log(q_upper),
            length.out = as.integer(config$binary_prior_q_grid_size %||% 9L))),
    config$pm,
    eb_q0,
    pmin(pmax(q_naive, q_lower), q_upper)
  )))
  pc_ped <- v22_binary_proband_only_pedigree(ped)

  scores <- lapply(q_grid, function(q) {
    obs <- v22_sample_binary_prior_genotypes(ped, q, config)
    pc <- v22_sample_binary_prior_genotypes(pc_ped, q, config)
    c(q = q,
      delta_A = obs$expected_A - pc$expected_A,
      eb_score_A = obs$expected_A - pc$expected_A + eb_n0 * (eb_q0 - q),
      score_q = obs$score_q_mean - pc$score_q_mean,
      obs_A = obs$expected_A,
      pc_A = pc$expected_A,
      obs_info_q = obs$info_q,
      pc_info_q = pc$info_q,
      founder_count = obs$founder_count)
  })
  score_mat <- do.call(rbind, scores)
  best <- which.min(abs(score_mat[, "eb_score_A"]))
  qhat <- as.numeric(score_mat[best, "q"])
  prior_info_q <- eb_n0 * eb_q0 / qhat^2 + eb_n0 * (1 - eb_q0) / (1 - qhat)^2
  info_q <- max(as.numeric(score_mat[best, "obs_info_q"] - score_mat[best, "pc_info_q"]) +
                  prior_info_q, 1e-8)
  var_q <- 1 / info_q
  list(
    q = qhat,
    vcov = var_q,
    convergence = TRUE,
    failure_reason = NA_character_,
    diagnostics = list(
      n_observed_nonproband = n_obs,
      q_naive_nonproband = q_naive,
      q_grid = q_grid,
      q_score_delta_A = score_mat[, "delta_A"],
      q_score_eb_A = score_mat[, "eb_score_A"],
      q_score = score_mat[, "score_q"],
      obs_info_q_grid = score_mat[, "obs_info_q"],
      pc_info_q_grid = score_mat[, "pc_info_q"],
      pc_expected_A_grid = score_mat[, "pc_A"],
      expected_founder_A_observed = score_mat[best, "obs_A"],
      expected_founder_A_proband_only = score_mat[best, "pc_A"],
      eb = isTRUE(eb),
      eb_q0 = eb_q0,
      eb_n0 = eb_n0,
      prior_info_q = prior_info_q,
      observed_info_q = info_q
    )
  )
}

v22_interpolate_binary_pc_A_to_q <- function(A, q_grid, pc_A_grid, config = v22_default_config(),
                                            eb = FALSE, eb_q0 = config$binary_eb_q0 %||% config$pm,
                                            eb_n0 = config$binary_eb_n0 %||% 0) {
  ok <- is.finite(q_grid) & is.finite(pc_A_grid)
  q_grid <- as.numeric(q_grid[ok])
  pc_A_grid <- as.numeric(pc_A_grid[ok])
  if (length(q_grid) < 2L || length(unique(pc_A_grid)) < 2L) {
    return(pmin(pmax(A / max(2, 2 * max(pc_A_grid, 1)), config$binary_q_lower), config$binary_q_upper))
  }
  score_grid <- A - pc_A_grid
  if (isTRUE(eb)) score_grid <- score_grid + max(as.numeric(eb_n0), 0) * (eb_q0 - q_grid)
  best <- which.min(abs(score_grid))
  if (is.finite(score_grid[best]) && abs(score_grid[best]) < 1e-8) return(q_grid[best])
  if (isTRUE(eb)) return(q_grid[best])
  ord <- order(pc_A_grid, q_grid)
  pc_A_grid <- pc_A_grid[ord]
  q_grid <- q_grid[ord]
  keep <- !duplicated(pc_A_grid)
  stats::approx(pc_A_grid[keep], q_grid[keep], xout = A,
                rule = 2, ties = "ordered")$y
}

v22_fit_binary_hwe_prior_from_genotype_draw <- function(G, ped_dat, reference_fit,
                                                       config = v22_default_config(),
                                                       eb = FALSE) {
  blocks <- v22_family_blocks(ped_dat)
  cache <- v22_build_pedigree_cache(ped_dat)
  founder_count <- max(v22_binary_founder_count(cache), 1)
  A <- v22_binary_founder_A_count(G, cache, blocks)
  diag <- reference_fit$diagnostics %||% list()
  q_grid <- diag$q_grid %||% numeric(0)
  pc_A_grid <- diag$pc_expected_A_grid %||% numeric(0)
  eb <- isTRUE(eb) || isTRUE(diag$eb)
  eb_q0 <- diag$eb_q0 %||% (config$binary_eb_q0 %||% config$pm)
  eb_n0 <- diag$eb_n0 %||% if (isTRUE(eb)) (config$binary_eb_n0 %||% 0) else 0
  qhat <- if (length(q_grid) && length(pc_A_grid)) {
    v22_interpolate_binary_pc_A_to_q(A, q_grid, pc_A_grid, config,
                                    eb = eb, eb_q0 = eb_q0, eb_n0 = eb_n0)
  } else {
    A / (2 * founder_count)
  }
  qhat <- pmin(pmax(qhat, config$binary_q_lower %||% 1e-4), config$binary_q_upper %||% 0.25)

  pc_info <- if (length(q_grid) && length(diag$pc_info_q_grid %||% numeric(0))) {
    stats::approx(q_grid, diag$pc_info_q_grid, xout = qhat, rule = 2, ties = "ordered")$y
  } else 0
  complete_info <- A / qhat^2 + (2 * founder_count - A) / (1 - qhat)^2
  prior_info_q <- if (isTRUE(eb)) {
    max(as.numeric(eb_n0), 0) * eb_q0 / qhat^2 +
      max(as.numeric(eb_n0), 0) * (1 - eb_q0) / (1 - qhat)^2
  } else 0
  info_q <- max(complete_info - pc_info + prior_info_q, 1e-8)
  list(prior_version = "B-R", eta = c(q = qhat),
       vcov = matrix(1 / info_q, 1, 1, dimnames = list("q", "q")),
       q = qhat, known = FALSE, convergence = TRUE,
       failure_reason = NA_character_,
       diagnostics = list(founder_A_count = A,
                          founder_count = founder_count,
                          complete_info_q = complete_info,
                          pc_info_q = pc_info,
                          eb = isTRUE(eb),
                          eb_q0 = eb_q0,
                          eb_n0 = eb_n0,
                          prior_info_q = prior_info_q,
                          observed_info_q = info_q,
                          source = "latent_genotype_draw"))
}

v22_fit_binary_hwe_prior <- function(dat, prior_version = c("B-O", "B-R"),
                                    config = v22_default_config(),
                                    pedigree_dat = NULL,
                                    eb = FALSE) {
  prior_version <- match.arg(prior_version)
  if (prior_version == "B-O") {
    return(list(prior_version = "B-O", eta = numeric(0), vcov = matrix(0, 0, 0),
                q = config$pm, known = TRUE, convergence = TRUE,
                failure_reason = NA_character_))
  }
  q_fit <- tryCatch(
    v22_estimate_binary_q_carrier_conditioned(dat, pedigree_dat, config, eb = eb),
    error = function(e) list(q = config$pm, vcov = 1e6, convergence = FALSE,
                             failure_reason = conditionMessage(e),
                             diagnostics = list())
  )
  if (!isTRUE(q_fit$convergence)) {
    return(list(prior_version = "B-R", eta = c(q = config$pm),
                vcov = matrix(q_fit$vcov, 1, 1, dimnames = list("q", "q")),
                q = config$pm, known = FALSE, convergence = FALSE,
                failure_reason = q_fit$failure_reason,
                diagnostics = q_fit$diagnostics))
  }
  qhat <- q_fit$q
  list(prior_version = "B-R", eta = c(q = qhat),
       vcov = matrix(q_fit$vcov, 1, 1, dimnames = list("q", "q")),
       q = qhat, known = FALSE, convergence = TRUE,
       failure_reason = NA_character_,
       diagnostics = q_fit$diagnostics)
}

v22_prepare_binary_pedigree <- function(dat, pedigree_dat = NULL) {
  if (is.null(pedigree_dat)) {
    ped <- dat
    match_idx <- seq_len(nrow(dat))
  } else {
    ped <- pedigree_dat
    if (!"mgene" %in% names(ped)) ped$mgene <- NA_real_
    ped$mgene <- NA_real_
    match_idx <- match(as.character(dat$indID), as.character(ped$indID))
    if (anyNA(match_idx)) stop("Analysis records are not all present in pedigree_dat.")
    ped$mgene[match_idx] <- dat$mgene
  }
  ped$mgene[as.numeric(ped$proband) == 1] <- 1L
  list(pedigree_dat = ped, analysis_match = match_idx)
}

v22_make_binary_update_arrays <- function(ped_dat, analysis_dat, analysis_match,
                                         U_analysis, H_analysis, omega) {
  n <- nrow(ped_dat)
  th <- v22_theta_from_omega(omega)
  out <- list(
    time = rep(0, n),
    t0 = rep(0, n),
    status = rep(0, n),
    newx = rep(0, n),
    U = rep(0, n),
    H = rep(0, n),
    base_eta = rep(0, n),
    is_analysis = rep(FALSE, n)
  )
  out$time[analysis_match] <- analysis_dat$time
  out$t0[analysis_match] <- analysis_dat$t0
  out$status[analysis_match] <- analysis_dat$status
  out$newx[analysis_match] <- analysis_dat$newx
  out$U[analysis_match] <- U_analysis
  out$H[analysis_match] <- H_analysis
  out$base_eta[analysis_match] <- th$beta_c * analysis_dat$newx + U_analysis
  out$is_analysis[analysis_match] <- TRUE
  out
}

v22_local_genotype_weights <- function(states, loc, Gfam, cache, q, th, arrays, fam_global_idx) {
  w <- numeric(length(states))
  for (s in seq_along(states)) {
    g <- states[s]
    if (is.na(cache$father[loc]) || is.na(cache$mother[loc])) {
      prior <- v22_founder_prob(g, q)
    } else {
      prior <- v22_mendel_prob(g, Gfam[cache$mother[loc]], Gfam[cache$father[loc]])
    }
    if (prior > 0 && length(cache$children[[loc]])) {
      for (ch in cache$children[[loc]]) {
        if (is.na(cache$mother[ch]) || is.na(cache$father[ch])) next
        mg <- if (cache$mother[ch] == loc) g else Gfam[cache$mother[ch]]
        fg <- if (cache$father[ch] == loc) g else Gfam[cache$father[ch]]
        prior <- prior * v22_mendel_prob(Gfam[ch], mg, fg)
      }
    }
    z <- as.integer(g >= 1L)
    global <- fam_global_idx[loc]
    disease <- arrays$status[global] * th$beta_b * z -
      arrays$H[global] * v22_safe_exp(arrays$base_eta[global] + th$beta_b * z)
    w[s] <- prior * v22_safe_exp(disease)
  }
  if (!any(is.finite(w)) || sum(w) <= 0) w <- rep(1, length(states))
  w / sum(w)
}

v22_draw_binary_genotype_sweep <- function(ped_dat, G, q, omega, arrays, config = v22_default_config(),
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
  for (b in seq_along(blocks)) {
    idx <- blocks[[b]]
    Gfam <- Gnew[idx]
    for (loc in seq_along(idx)) {
      states <- states_list[[idx[loc]]]
      if (length(states) == 1L) {
        Gfam[loc] <- states[1]
      } else {
        w <- v22_local_genotype_weights(states, loc, Gfam, cache[[b]], q, th, arrays, idx)
        Gfam[loc] <- sample(states, 1L, prob = w)
      }
    }
    Gnew[idx] <- Gfam
  }
  Gnew
}

v22_draw_binary_chain <- function(dat, K, omega, prior_fit, pedigree_dat = NULL,
                                 M = 10L, burnin = 50L, thin = 10L,
                                 config = v22_default_config(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- v22_align_K(K, dat)
  ped <- v22_prepare_binary_pedigree(dat, pedigree_dat)
  ped_dat <- ped$pedigree_dat
  analysis_match <- ped$analysis_match
  q <- prior_fit$q
  G <- v22_initialize_binary_genotypes(ped_dat, q)
  blocks_info <- v22_precompute_K_blocks(K, dat)
  ped_blocks <- v22_family_blocks(ped_dat)
  ped_cache <- v22_build_pedigree_cache(ped_dat)
  zobs <- as.numeric(ped_dat$mgene)
  is_prob <- as.numeric(ped_dat$proband) == 1
  states_list <- lapply(seq_len(nrow(ped_dat)), function(i) {
    v22_admissible_genotypes(zobs[i], is_prob[i])
  })
  H_analysis <- v22_H0_diff(dat$time, dat$t0, omega, config$agemin)
  th <- v22_theta_from_omega(omega)
  completed <- vector("list", M)
  draws <- vector("list", M)
  save_i <- 0L
  total <- burnin + M * thin
  for (iter in seq_len(total)) {
    z_analysis <- as.integer(G[analysis_match] >= 1L)
    dat_for_u <- dat
    dat_for_u$mgene <- z_analysis
    U <- v22_draw_u_vg(dat_for_u, K, omega, dat_for_u$newx, config,
                      blocks_info = blocks_info)$u
    arrays <- v22_make_binary_update_arrays(ped_dat, dat, analysis_match, U,
                                           H_analysis, omega)
    for (sweep in seq_len(config$binary_genotype_sweeps)) {
      G <- v22_draw_binary_genotype_sweep(ped_dat, G, q, omega, arrays, config,
                                         blocks = ped_blocks, cache = ped_cache,
                                         states_list = states_list, th = th)
    }
    if (iter > burnin && ((iter - burnin) %% thin == 0L)) {
      save_i <- save_i + 1L
      z_analysis <- as.integer(G[analysis_match] >= 1L)
      di <- dat
      di$mgene <- z_analysis
      di$mgene[as.numeric(di$proband) == 1] <- 1L
      completed[[save_i]] <- di
      draws[[save_i]] <- list(mgene = z_analysis, G = G, U = U, q = q)
    }
  }
  z_draw_mat <- if (length(draws)) {
    do.call(cbind, lapply(draws, function(x) as.integer(x$mgene)))
  } else {
    matrix(integer(0), nrow = nrow(dat), ncol = 0)
  }
  missing_analysis <- is.na(dat$mgene)
  posterior_carrier_prob <- rep(NA_real_, nrow(dat))
  if (ncol(z_draw_mat)) {
    posterior_carrier_prob[missing_analysis] <- rowMeans(z_draw_mat[missing_analysis, , drop = FALSE])
  }
  p_mis <- posterior_carrier_prob[missing_analysis]
  p_mis_finite <- p_mis[is.finite(p_mis)]
  list(completed = completed, draws = draws,
       diagnostics = list(n_missing = sum(is.na(dat$mgene)),
                          missing_proband_carrier = sum(is.na(dat$mgene) & as.numeric(dat$proband) == 1),
                          posterior_carrier_prob_missing = p_mis,
                          posterior_carrier_prob_mean = if (length(p_mis_finite)) mean(p_mis_finite) else NA_real_,
                          posterior_carrier_prob_max = if (length(p_mis_finite)) max(p_mis_finite) else NA_real_))
}

v22_draw_binary_prior_posterior <- function(prior_fit, config = v22_default_config()) {
  if (prior_fit$prior_version == "B-O") return(prior_fit)
  if (!isTRUE(prior_fit$convergence)) {
    stop("Cannot draw B-R allele frequency from a failed completed-data prior fit.")
  }
  q <- as.numeric(prior_fit$q %||% prior_fit$eta["q"])
  q <- pmin(pmax(q, config$binary_q_lower %||% 1e-4), config$binary_q_upper %||% 0.25)
  Vq <- as.numeric(prior_fit$vcov["q", "q"])
  if (!is.finite(Vq) || Vq <= 0) stop("B-R q covariance is not positive.")
  phi <- qlogis(q)
  Vphi <- Vq / max(q * (1 - q), 1e-8)^2
  q_draw <- v22_inv_logit(stats::rnorm(1L, phi, sqrt(max(Vphi, 1e-12))))
  q_draw <- pmin(pmax(q_draw, config$binary_q_lower %||% 1e-4), config$binary_q_upper %||% 0.25)
  out <- prior_fit
  out$q <- q_draw
  out$eta <- c(q = q_draw)
  out$posterior_draw <- TRUE
  out
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
  K <- v22_align_K(K, dat)
  ped <- v22_prepare_binary_pedigree(dat, pedigree_dat)
  ped_dat <- ped$pedigree_dat
  analysis_match <- ped$analysis_match
  use_beta_penalty <- identical(prior_version, "B-R")
  base_prior <- v22_fit_binary_hwe_prior(dat, prior_version, config,
                                         pedigree_dat = pedigree_dat,
                                         eb = use_beta_penalty)
  if (prior_version == "B-R" && !isTRUE(base_prior$convergence)) {
    v22_stop_pdmi_diagnostic(
      base_prior$failure_reason %||% "Initial B-R prior fit failed.",
      list(stage = c(debug_context, list(phase = "initial_prior_fit")),
           prior_fit = base_prior)
    )
  }
  coef_names <- c(v22_omega_names(), if (prior_version == "B-R") "q" else character(0))
  coefIter <- array(NA_real_, dim = c(M, length(coef_names), numit),
                    dimnames = list(imputation = seq_len(M), parameter = coef_names,
                                    iteration = seq_len(numit)))
  completed <- vector("list", M)
  draws <- vector("list", M)
  blocks_info <- v22_precompute_K_blocks(K, dat)
  ped_blocks <- v22_family_blocks(ped_dat)
  ped_cache <- v22_build_pedigree_cache(ped_dat)
  zobs <- as.numeric(ped_dat$mgene)
  is_prob <- as.numeric(ped_dat$proband) == 1
  states_list <- lapply(seq_len(nrow(ped_dat)), function(i) {
    v22_admissible_genotypes(zobs[i], is_prob[i])
  })
  if (!"t0" %in% names(dat)) dat$t0 <- 0
  trace <- list()
  add_trace <- function(entry) {
    if (v22_pdmi_debug_enabled(config)) {
      trace[[length(trace) + 1L]] <<- entry
      if (length(trace) > 25L) trace <<- trace[(length(trace) - 24L):length(trace)]
    }
    invisible(NULL)
  }
  for (m in seq_len(M)) {
    prior_current <- base_prior
    G <- v22_initialize_binary_genotypes(ped_dat, prior_current$q)
    chain_init_omega <- init_omega
    iter_draws <- vector("list", numit)
    for (iter in seq_len(numit)) {
      z_analysis <- as.integer(G[analysis_match] >= 1L)
      current <- dat
      current$mgene <- z_analysis
      current$mgene[as.numeric(current$proband) == 1] <- 1L
      stage <- c(debug_context, list(missing_type = "binary",
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
      prior_fit <- if (prior_version == "B-R") {
        v22_fit_binary_hwe_prior_from_genotype_draw(G, ped_dat, base_prior,
                                                    config, eb = TRUE)
      } else {
        base_prior
      }
      prior_draw <- tryCatch(v22_draw_binary_prior_posterior(prior_fit, config),
                             error = function(e) e)
      if (inherits(prior_draw, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(prior_draw),
          list(stage = c(stage, list(phase = "prior_draw")),
               trace_tail = trace,
               prior_fit = prior_fit)
        )
      }
      H_analysis <- v22_H0_diff(dat$time, dat$t0, omega_draw, config$agemin)
      U_draw <- tryCatch(v22_draw_u_vg(current, K, omega_draw, current$newx, config,
                                       blocks_info = blocks_info),
                         error = function(e) e)
      if (inherits(U_draw, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(U_draw),
          list(stage = c(stage, list(phase = "binary_u_draw")),
               trace_tail = trace,
               current_omega = disease_fit$omega,
               omega_draw = omega_draw,
               root_cause = v22_classify_pdmi_failure(conditionMessage(U_draw),
                                                       list(omega_draw = omega_draw)))
        )
      }
      U <- U_draw$u
      arrays <- v22_make_binary_update_arrays(ped_dat, dat, analysis_match, U,
                                              H_analysis, omega_draw)
      th <- v22_theta_from_omega(omega_draw)
      genotype_update <- tryCatch({
        for (sweep in seq_len(config$binary_genotype_sweeps)) {
          G <- v22_draw_binary_genotype_sweep(ped_dat, G, prior_draw$q, omega_draw, arrays,
                                              config, blocks = ped_blocks, cache = ped_cache,
                                              states_list = states_list, th = th)
        }
        G
      }, error = function(e) e)
      if (inherits(genotype_update, "error")) {
        v22_stop_pdmi_diagnostic(
          conditionMessage(genotype_update),
          list(stage = c(stage, list(phase = "binary_genotype_update")),
               trace_tail = trace,
               current_omega = disease_fit$omega,
               omega_draw = omega_draw,
               q_draw = prior_draw$q)
        )
      }
      G <- genotype_update
      coefIter[m, v22_omega_names(), iter] <- omega_draw[v22_omega_names()]
      if (prior_version == "B-R") coefIter[m, "q", iter] <- prior_draw$q
      chain_init_omega <- disease_fit$omega
      iter_draws[[iter]] <- list(omega = omega_draw, q = prior_draw$q)
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
         beta_prior_q0 = config$binary_beta_prior_q0 %||% config$pm,
         beta_prior_n0 = config$binary_beta_prior_n0 %||% 0,
         trace_tail = trace
       ))
}

v22_binary_systematic_counts <- function(p, M) {
  p <- as.numeric(p)
  p[!is.finite(p) | p < 0] <- 0
  if (!(sum(p) > 0)) p <- rep(1 / length(p), length(p)) else p <- p / sum(p)
  n0 <- floor(M * p)
  R <- M - sum(n0)
  if (R <= 0L) return(n0)
  residual <- M * p - n0
  if (!(sum(residual) > 0)) {
    n0[seq_len(R)] <- n0[seq_len(R)] + 1L
    return(n0)
  }
  r <- residual / sum(residual)
  u <- stats::runif(1, 0, 1 / R)
  grid <- u + (seq_len(R) - 1) / R
  cs <- cumsum(r)
  add <- tabulate(findInterval(grid, c(0, cs), rightmost.closed = TRUE), nbins = length(p))
  n0 + add
}

v22_binary_config_key <- function(z) paste(as.integer(z), collapse = "")

v22_binary_config_values <- function(key) {
  if (!nzchar(key)) return(integer(0))
  as.integer(strsplit(key, "", fixed = TRUE)[[1]])
}

v22_calibrate_binary_imputations <- function(dat, ped_dat, candidate_imp,
                                            M = 10L, config = v22_default_config(),
                                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M <- max(1L, as.integer(M))
  S <- length(candidate_imp$draws)
  if (!S) stop("No candidate binary draws available for calibration.")
  completed <- replicate(M, dat, simplify = FALSE)
  G_completed <- replicate(M, integer(nrow(ped_dat)), simplify = FALSE)
  fams <- unique(dat$famID)
  config_counts <- integer(0)

  for (fam in fams) {
    idx <- which(dat$famID == fam)
    ped_idx <- which(ped_dat$famID == fam)
    miss_idx <- idx[is.na(dat$mgene[idx])]
    if (!length(miss_idx)) {
      chosen <- sample.int(S, M, replace = TRUE)
      for (m in seq_len(M)) {
        G_completed[[m]][ped_idx] <- candidate_imp$draws[[chosen[m]]]$G[ped_idx]
      }
      next
    }
    keys <- vapply(candidate_imp$draws, function(draw) {
      v22_binary_config_key(draw$mgene[miss_idx])
    }, character(1))
    tab <- table(keys)
    uniq <- names(tab)
    counts <- v22_binary_systematic_counts(as.numeric(tab) / sum(tab), M)
    assigned <- rep(uniq, counts)
    if (length(assigned) != M) {
      assigned <- sample(uniq, M, replace = TRUE, prob = as.numeric(tab))
    } else {
      assigned <- sample(assigned, M, replace = FALSE)
    }
    config_counts <- c(config_counts, length(uniq))
    for (m in seq_len(M)) {
      vals <- v22_binary_config_values(assigned[m])
      completed[[m]]$mgene[miss_idx] <- vals
      matched_draws <- which(keys == assigned[m])
      chosen <- sample(matched_draws, 1L)
      G_completed[[m]][ped_idx] <- candidate_imp$draws[[chosen]]$G[ped_idx]
    }
  }
  for (m in seq_len(M)) {
    completed[[m]]$mgene[as.numeric(completed[[m]]$proband) == 1] <- 1L
  }
  z_draw_mat <- do.call(cbind, lapply(completed, function(x) as.integer(x$mgene)))
  missing_analysis <- is.na(dat$mgene)
  p_mis <- if (any(missing_analysis)) rowMeans(z_draw_mat[missing_analysis, , drop = FALSE]) else numeric(0)
  list(
    completed = completed,
    draws = lapply(seq_len(M), function(m) list(
      mgene = as.integer(completed[[m]]$mgene),
      G = G_completed[[m]],
      q = candidate_imp$draws[[1]]$q
    )),
    diagnostics = list(
      calibration = "residual_systematic",
      candidate_draws = S,
      n_missing = sum(is.na(dat$mgene)),
      missing_proband_carrier = sum(is.na(dat$mgene) & as.numeric(dat$proband) == 1),
      families_with_missing = sum(vapply(split(is.na(dat$mgene), dat$famID), any, logical(1))),
      mean_unique_configs_per_incomplete_family = if (length(config_counts)) mean(config_counts) else 0,
      posterior_carrier_prob_missing = p_mis,
      posterior_carrier_prob_mean = if (length(p_mis)) mean(p_mis) else NA_real_,
      posterior_carrier_prob_max = if (length(p_mis)) max(p_mis) else NA_real_
    )
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
  K <- v22_align_K(K, dat_with_missing)
  init_fit <- v22_fit_mean_completed_initial(dat_with_missing, K, "binary", config)
  omega <- if (isTRUE(init_fit$convergence)) init_fit$omega else v22_actual_omega(config$sigma_u2_grid[1], config)
  M <- as.integer(config$M_imp_pdmi %||% config$M_imp_cong)
  imp <- tryCatch(v22_draw_binary_pdmi(dat_with_missing, K, prior_version, pedigree_dat,
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
                        context = list(missing_type = "binary",
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
                       m_success = pool$M_success)
  )
}
