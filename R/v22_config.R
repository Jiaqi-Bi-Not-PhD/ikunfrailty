## ============================================================
## V2.2 frailtypack-congenial simulation configuration
## Technical-note source: Simulation design, scenario grid, and
## reported scale omega = (log rho, log lambda, beta_b, beta_c, sigma_u^2).
## ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

v22_omega_names <- function() {
  c("log.rho", "log.lambda", "beta_b", "beta_c", "sigma_u2")
}

v22_default_config <- function() {
  list(
    B_sim = 200L,
    M_imp_pdmi = 20L,
    M_imp_cong = 20L,
    M_imp_smcfcs = 20L,
    n_families = 498L,
    sigma_u2_grid = c(0.5, 0.2),
    missing_rates = c(0.20, 0.50, 0.80),
    pm = 0.02,
    omega_base = c(log.rho = 0.804, log.lambda = 4.71,
                   beta_b = 2.2, beta_c = 1.0),
    design = "pop+",
    variation = "kinship",
    base.dist = "Weibull",
    frailty.dist = "lognormal",
    interaction = FALSE,
    probandage = c(45, 2),
    agemin = 0,
    agemax = 100,
    selection_max_attempts_per_family = 5000L,
    align_female_only = TRUE,
    analysis_female_only = TRUE,
    prs_sigma2 = 0.1,
    penetrance_ages = c(40, 50, 60, 70, 80),
    penetrance_prs = c(-0.5, 0, 0.5),
    penetrance_gene = c(0, 1),
    penetrance_k0 = 1,
    frailtypack_maxit = 35L,
    gh_order = 20L,
    mcmc_burnin = 50L,
    mcmc_thin = 10L,
    vg_maxit = 50L,
    vg_tol = 1e-5,
    continuous_slice_step = 0.7,
    continuous_slice_m = 80L,
    pdmi_numit = 10L,
    pdmi_prior_version = "exact_slice_weak_proper_v1",
    theta_prior_sd_log_baseline = 10,
    theta_prior_sd_beta = 10,
    theta_prior_tau_shape = 1.1,
    theta_prior_tau_rate = 0.2,
    theta_tau_lower = 1e-6,
    theta_tau_upper = 5,
    theta_slice_widths = c(log.rho = 0.08, log.lambda = 0.08,
                           beta_b = 0.12, beta_c = 0.12,
                           sigma_u2 = 0.05),
    theta_slice_sweeps = 1L,
    theta_slice_m = 80L,
    ess_max_shrink = 200L,
    cr_gamma_prior_sd = 10,
    cr_log_var_prior_sd = 3,
    cr_log_var_slice_width = 0.5,
    cr_log_var_slice_m = 60L,
    frailtypack_sigma_lower = 1e-8,
    binary_genotype_sweeps = 8L,
    binary_prior_q_grid_size = 9L,
    binary_prior_gibbs_burnin = 30L,
    binary_prior_gibbs_draws = 20L,
    binary_prior_gibbs_thin = 2L,
    binary_q_lower = 1e-4,
    binary_q_upper = 0.25,
    binary_beta_prior_q0 = 0.02,
    binary_beta_prior_n0 = 50,
    binary_eb_q0 = 0.02,
    binary_eb_n0 = 50,
    pdmi_debug = tolower(Sys.getenv("SIM_PDMI_DEBUG", "0")) %in% c("1", "true", "yes", "y"),
    mcem_iter = 2L,
    mcem_draws = 3L,
    mcem_burnin = 20L,
    skip_existing_results = TRUE,
    skip_existing_benchmarks = TRUE,
    run_label = Sys.getenv("SIM_RUN_LABEL", paste0(Sys.info()[["user"]], "_", format(Sys.time(), "%Y%m%d_%H%M%S"))),
    results_root = Sys.getenv("SIM_RESULTS_ROOT", file.path("Results", "raw"))
  )
}

v22_actual_omega <- function(sigma_u2, config = v22_default_config()) {
  c(config$omega_base, sigma_u2 = sigma_u2)[v22_omega_names()]
}

v22_actual_value_legacy_names <- function(sigma_u2, config = v22_default_config()) {
  c(log.shape = unname(config$omega_base["log.rho"]),
    log.scale = unname(config$omega_base["log.lambda"]),
    beta_mgene = unname(config$omega_base["beta_b"]),
    beta_PRS = unname(config$omega_base["beta_c"]),
    sigma2 = sigma_u2)
}

v22_theta_from_omega <- function(omega) {
  omega <- omega[v22_omega_names()]
  list(
    log_rho = unname(omega["log.rho"]),
    rho = exp(unname(omega["log.rho"])),
    log_lambda = unname(omega["log.lambda"]),
    lambda = exp(unname(omega["log.lambda"])),
    beta_b = unname(omega["beta_b"]),
    beta_c = unname(omega["beta_c"]),
    sigma_u2 = unname(omega["sigma_u2"])
  )
}

v22_log_kappa_from_omega <- function(omega) {
  th <- v22_theta_from_omega(omega)
  -th$rho * th$log_lambda
}

v22_method_tag <- function(method, prior_version = NA_character_) {
  method <- tolower(method)
  prior_version <- tolower(prior_version %||% "")
  if (identical(method, "mi-smcfcs")) return("mi_smcfcs")
  if (identical(method, "c-o-pdmi")) return("c_o_pdmi")
  if (identical(method, "c-r-pdmi")) return("c_r_pdmi")
  if (identical(method, "b-o-pdmi")) return("b_o_pdmi")
  if (identical(method, "b-r-pdmi")) return("b_r_pdmi")
  if (identical(method, "cca")) return("cca")
  if (identical(method, "j-smcfcs")) return("j_smcfcs")
  if (identical(method, "j-o-pdmi")) return("j_o_pdmi")
  if (identical(method, "j-r-pdmi")) return("j_r_pdmi")
  gsub("[^a-z0-9]+", "_", paste(method, prior_version))
}

v22_parse_numeric_env <- function(name, default) {
  value <- Sys.getenv(name, paste(default, collapse = ","))
  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  pieces <- pieces[nzchar(pieces)]
  out <- suppressWarnings(as.numeric(pieces))
  if (!length(out) || anyNA(out)) {
    stop("Environment variable ", name, " must contain numeric value(s), got: ", value)
  }
  out
}

v22_parse_single_numeric_env <- function(name, default, label = name) {
  values <- v22_parse_numeric_env(name, default)
  if (length(values) != 1L) {
    stop(label, " must contain exactly one value per job. Got ",
         name, "=", paste(values, collapse = ","),
         ". Submit separate jobs for repeated grid values.")
  }
  values[[1]]
}

v22_number_tag <- function(prefix, x) {
  paste0(prefix, "_", gsub("\\.", "p", as.character(x)))
}

v22_rate_tag <- function(x) {
  if (is.na(x)) return("nomiss")
  sprintf("miss%02d", round(100 * x))
}

v22_clean_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

v22_find_code_root <- function(starts = c(getwd())) {
  for (start in starts) {
    if (!nzchar(start)) next
    cur <- normalizePath(start, mustWork = FALSE)
    if (!dir.exists(cur)) cur <- dirname(cur)
    repeat {
      if (dir.exists(file.path(cur, "Shared")) &&
          dir.exists(file.path(cur, "Continuous Only")) &&
          dir.exists(file.path(cur, "Binary Only"))) {
        return(normalizePath(cur))
      }
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
  }
  stop("Could not locate the V2.2 code root.")
}

v22_get_job_cores <- function(default = 1L) {
  env_names <- c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE")
  for (name in env_names) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      first <- strsplit(value, ",", fixed = TRUE)[[1]][1]
      n <- suppressWarnings(as.integer(sub("\\(.*$", "", first)))
      if (!is.na(n) && n > 0L) return(n)
    }
  }
  n <- suppressWarnings(parallel::detectCores())
  if (!is.na(n) && n > 0L) return(n)
  as.integer(default)
}

v22_seed_streams <- function(replicate_id, missing_type = "none", target_missing_rate = NA_real_,
                            method = "full-data", seed_base = 700000L) {
  type_offset <- switch(missing_type,
                        none = 0L,
                        continuous = 10000L,
                        binary = 20000L,
                        30000L)
  rate_offset <- if (is.na(target_missing_rate)) 0L else as.integer(round(1000 * target_missing_rate))
  method_offset <- sum(utf8ToInt(method)) %% 10000L
  list(
    complete_seed = seed_base + replicate_id,
    missing_mask_seed = seed_base + 100000L + type_offset + rate_offset + replicate_id,
    method_seed = seed_base + 200000L + type_offset + rate_offset + method_offset + replicate_id
  )
}

v22_result_metadata <- function(replicate_id, missing_type, sigma_u2, target_missing_rate,
                               method, prior_version, M_imp, seeds, config) {
  list(
    replicate_id = replicate_id,
    missing_type = missing_type,
    sigma_u2 = sigma_u2,
    target_missing_rate = target_missing_rate,
    method = method,
    prior_version = prior_version,
    pdmi_prior_version = config$pdmi_prior_version %||% "exact_slice_weak_proper_v1",
    M_imp = M_imp,
    complete_data_seed = seeds$complete_seed,
    missing_mask_seed = seeds$missing_mask_seed,
    method_seed = seeds$method_seed,
    run_label = config$run_label,
    run_owner = Sys.info()[["user"]],
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}
