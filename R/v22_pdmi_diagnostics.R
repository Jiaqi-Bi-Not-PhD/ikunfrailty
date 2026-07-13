## ============================================================
## PDMI finite-value and frailtypack failure diagnostics.
## These helpers are diagnostic-only: they never drop rows, clip
## estimates, or change the target analysis dataset.
## ============================================================

v22_pdmi_debug_enabled <- function(config = v22_default_config()) {
  isTRUE(config$pdmi_debug %||% FALSE) ||
    tolower(Sys.getenv("SIM_PDMI_DEBUG", "0")) %in% c("1", "true", "yes", "y")
}

v22_numeric_summary <- function(x, probs = c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1)) {
  x <- suppressWarnings(as.numeric(x))
  finite <- x[is.finite(x)]
  out <- setNames(rep(NA_real_, length(probs)), paste0("q", gsub("\\.", "_", probs)))
  if (length(finite)) {
    out <- stats::quantile(finite, probs = probs, na.rm = TRUE, names = FALSE, type = 7)
    names(out) <- paste0("q", gsub("\\.", "_", probs))
  }
  c(n = length(x),
    n_finite = sum(is.finite(x)),
    n_na = sum(is.na(x)),
    n_nan = sum(is.nan(x)),
    n_pos_inf = sum(is.infinite(x) & x > 0),
    n_neg_inf = sum(is.infinite(x) & x < 0),
    out)
}

v22_column_finite_diagnostics <- function(dat, columns, finite_columns = columns) {
  rows <- lapply(columns, function(col) {
    present <- col %in% names(dat)
    if (!present) {
      return(data.frame(column = col, present = FALSE, n = NA_integer_,
                        n_finite = NA_integer_, n_na = NA_integer_,
                        n_nan = NA_integer_, n_pos_inf = NA_integer_,
                        n_neg_inf = NA_integer_, n_nonfinite = NA_integer_,
                        q0 = NA_real_, q0_01 = NA_real_, q0_05 = NA_real_,
                        q0_5 = NA_real_, q0_95 = NA_real_, q0_99 = NA_real_,
                        q1 = NA_real_, stringsAsFactors = FALSE))
    }
    x <- dat[[col]]
    if (col %in% finite_columns) {
      s <- v22_numeric_summary(x)
      n_nonfinite <- as.integer(s["n"] - s["n_finite"])
      return(data.frame(column = col, present = TRUE,
                        n = as.integer(s["n"]),
                        n_finite = as.integer(s["n_finite"]),
                        n_na = as.integer(s["n_na"]),
                        n_nan = as.integer(s["n_nan"]),
                        n_pos_inf = as.integer(s["n_pos_inf"]),
                        n_neg_inf = as.integer(s["n_neg_inf"]),
                        n_nonfinite = n_nonfinite,
                        q0 = s["q0"], q0_01 = s["q0_01"], q0_05 = s["q0_05"],
                        q0_5 = s["q0_5"], q0_95 = s["q0_95"],
                        q0_99 = s["q0_99"], q1 = s["q1"],
                        stringsAsFactors = FALSE))
    }
    data.frame(column = col, present = TRUE,
               n = length(x), n_finite = NA_integer_,
               n_na = sum(is.na(x)), n_nan = NA_integer_,
               n_pos_inf = NA_integer_, n_neg_inf = NA_integer_,
               n_nonfinite = sum(is.na(x)),
               q0 = NA_real_, q0_01 = NA_real_, q0_05 = NA_real_,
               q0_5 = NA_real_, q0_95 = NA_real_, q0_99 = NA_real_,
               q1 = NA_real_, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

v22_carrier_event_diagnostics <- function(dat) {
  if (!all(c("mgene", "status") %in% names(dat))) return(data.frame())
  z <- suppressWarnings(as.integer(as.numeric(dat$mgene) >= 1L))
  d <- suppressWarnings(as.integer(as.numeric(dat$status) == 1L))
  ok <- is.finite(z) & is.finite(d)
  if (!any(ok)) return(data.frame())
  tab <- as.data.frame(table(carrier = z[ok], event = d[ok]), stringsAsFactors = FALSE)
  names(tab)[names(tab) == "Freq"] <- "n"
  tab$carrier <- as.integer(as.character(tab$carrier))
  tab$event <- as.integer(as.character(tab$event))
  tab
}

v22_K_alignment_diagnostics <- function(K, dat, kinship_cache = NULL) {
  out <- list(
    n_dat = nrow(dat),
    n_K = if (is.null(K)) NA_integer_ else nrow(K),
    k_square = !is.null(K) && is.matrix(as.matrix(K)) && nrow(K) == ncol(K),
    rownames_present = !is.null(rownames(K)),
    align_ok = FALSE,
    align_error = NA_character_,
    n_families = if ("famID" %in% names(dat)) length(unique(dat$famID)) else NA_integer_,
    first_bad_family = NA_character_,
    min_family_eigen = NA_real_,
    cache_validated = FALSE
  )
  if (v22_kinship_cache_matches(kinship_cache, dat)) {
    out$n_K <- nrow(kinship_cache$K)
    out$k_square <- TRUE
    out$rownames_present <- TRUE
    out$align_ok <- TRUE
    out$n_families <- kinship_cache$n_families %||% out$n_families
    out$cache_validated <- TRUE
    return(out)
  }
  aligned <- tryCatch(v22_align_K(K, dat), error = function(e) e)
  if (inherits(aligned, "error")) {
    out$align_error <- conditionMessage(aligned)
    return(out)
  }
  out$align_ok <- TRUE
  blocks <- tryCatch(v22_family_blocks(dat), error = function(e) list())
  minev <- Inf
  for (nm in names(blocks)) {
    idx <- blocks[[nm]]
    Ki <- aligned[idx, idx, drop = FALSE]
    ev <- tryCatch(eigen(0.5 * (Ki + t(Ki)), symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    if (!all(is.finite(ev))) {
      out$first_bad_family <- nm
      out$min_family_eigen <- NA_real_
      return(out)
    }
    minev <- min(minev, min(ev))
    if (min(ev) < -1e-6 && is.na(out$first_bad_family)) out$first_bad_family <- nm
  }
  out$min_family_eigen <- if (is.finite(minev)) minev else NA_real_
  out
}

v22_analysis_data_diagnostics <- function(dat, K = NULL, context = list(),
                                          kinship_cache = NULL) {
  required <- c("t0", "time", "status", "mgene", "newx", "famID",
                "proband", "currentage", "indID")
  finite_cols <- c("t0", "time", "status", "mgene", "newx", "proband", "currentage")
  cols <- v22_column_finite_diagnostics(dat, required, finite_cols)
  bad_cols <- cols[!cols$present | cols$n_nonfinite > 0, , drop = FALSE]
  list(
    context = context,
    n_rows = nrow(dat),
    column_diagnostics = cols,
    bad_columns = bad_cols$column,
    has_bad_input = nrow(bad_cols) > 0,
    has_na = any(cols$n_na %||% 0, na.rm = TRUE),
    has_inf = any((cols$n_pos_inf %||% 0) > 0 | (cols$n_neg_inf %||% 0) > 0, na.rm = TRUE),
    newx_summary = if ("newx" %in% names(dat)) v22_numeric_summary(dat$newx) else numeric(0),
    carrier_event = v22_carrier_event_diagnostics(dat),
    K = if (!is.null(K)) v22_K_alignment_diagnostics(K, dat, kinship_cache) else list()
  )
}

v22_pdmi_diagnostic_error <- function(message, diagnostics = list()) {
  structure(list(message = message, call = NULL, diagnostics = diagnostics),
            class = c("v22_pdmi_diagnostic_error", "error", "condition"))
}

v22_pdmi_stage_phase <- function(stage, phase) {
  stage <- stage[setdiff(names(stage), "phase")]
  c(stage, list(phase = phase))
}

v22_stop_pdmi_diagnostic <- function(message, diagnostics = list()) {
  stop(v22_pdmi_diagnostic_error(message, diagnostics))
}

v22_condition_diagnostics <- function(e) {
  if (inherits(e, "v22_pdmi_diagnostic_error")) return(e$diagnostics %||% list())
  list(error_message = conditionMessage(e), error_class = class(e)[1])
}

v22_omega_is_extreme <- function(omega) {
  if (is.null(omega) || any(!is.finite(omega[v22_omega_names()]))) return(TRUE)
  abs(omega["log.rho"]) > 10 ||
    abs(omega["log.lambda"]) > 10 ||
    abs(omega["beta_b"]) > 30 ||
    abs(omega["beta_c"]) > 30 ||
    omega["sigma_u2"] <= 1e-8 ||
    omega["sigma_u2"] > 50
}

v22_frailtypack_raw_diagnostics <- function(fit) {
  if (is.null(fit)) return(list())
  Vraw <- tryCatch(as.matrix(fit$varHtotal), error = function(e) NULL)
  b <- tryCatch(as.numeric(fit$b), error = function(e) numeric(0))
  list(
    istop = fit$istop %||% NA_integer_,
    shape_weib = tryCatch(as.numeric(fit$shape.weib), error = function(e) NA_real_),
    scale_weib = tryCatch(as.numeric(fit$scale.weib), error = function(e) NA_real_),
    coef = tryCatch(fit$coef, error = function(e) numeric(0)),
    sigma2 = tryCatch(as.numeric(fit$sigma2), error = function(e) NA_real_),
    b_length = length(b),
    b_all_finite = length(b) > 0L && all(is.finite(b)),
    b_nonfinite_count = sum(!is.finite(b)),
    varHtotal_dim = if (is.null(Vraw)) c(NA_integer_, NA_integer_) else dim(Vraw),
    varHtotal_all_finite = !is.null(Vraw) && all(is.finite(Vraw)),
    varHtotal_psd = !is.null(Vraw) && v22_is_psd(Vraw)
  )
}

v22_fit_status_table <- function(fits) {
  if (!length(fits)) return(data.frame())
  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    omega <- fit$omega %||% setNames(rep(NA_real_, 5), v22_omega_names())
    attempts <- fit$diagnostics$final_fit_init_attempts %||% list()
    attempt_labels <- vapply(attempts, function(x) x$label %||% NA_character_, character(1))
    attempt_ok <- vapply(attempts, function(x) isTRUE(x$convergence), logical(1))
    scenario_idx <- match("scenario_omega", attempt_labels)
    sampler_idx <- match("sampler_last_omega", attempt_labels)
    cbind(data.frame(
      m = i,
      convergence = isTRUE(fit$convergence),
      failure_reason = fit$failure_reason %||% NA_character_,
      selected_init_label = fit$diagnostics$selected_init_label %||% NA_character_,
      n_init_attempts = length(attempts),
      scenario_init_convergence = if (is.na(scenario_idx)) NA else attempt_ok[[scenario_idx]],
      sampler_last_init_convergence = if (is.na(sampler_idx)) NA else attempt_ok[[sampler_idx]],
      vcov_all_finite = !is.null(fit$vcov_omega) && all(is.finite(fit$vcov_omega)),
      stringsAsFactors = FALSE
    ), as.data.frame(as.list(omega[v22_omega_names()]), check.names = FALSE))
  })
  out <- do.call(rbind, rows)
  names(out) <- make.names(names(out), unique = TRUE)
  out
}

v22_classify_pdmi_failure <- function(failure_reason, diagnostics = list()) {
  reason <- failure_reason %||% ""
  finite_diag <- diagnostics$frailtypack$input_diagnostics %||%
    diagnostics$input_diagnostics %||% diagnostics$finite_diagnostics %||% list()
  cols <- finite_diag$column_diagnostics %||% data.frame()
  if (nrow(cols) && any(cols$n_nonfinite > 0, na.rm = TRUE)) return("input_nonfinite")
  if (grepl("Non-finite analysis data|Completed analysis data", reason, fixed = FALSE)) {
    return("input_nonfinite")
  }
  newx <- finite_diag$newx_summary %||% diagnostics$newx_summary %||% numeric(0)
  if (length(newx) && any(abs(newx[c("q0", "q1")] %||% 0) > 20, na.rm = TRUE)) {
    return("extreme_x")
  }
  if (grepl("alpha", reason, ignore.case = TRUE)) return("alpha_underflow")
  omega <- diagnostics$omega_draw %||% diagnostics$current_omega %||% numeric(0)
  if (length(omega) && v22_omega_is_extreme(omega)) return("omega_draw_extreme")
  raw <- diagnostics$frailtypack$raw %||% diagnostics$raw %||% list()
  sig <- raw$sigma2 %||% NA_real_
  if (length(sig) && (!is.finite(sig[1]) || sig[1] <= 1e-8)) return("frailty_boundary")
  if (grepl("PDMI pooling requires", reason, fixed = TRUE)) return("final_fit_only")
  carrier_event <- finite_diag$carrier_event %||% data.frame()
  if (nrow(carrier_event)) {
    events_by_carrier <- stats::aggregate(n ~ carrier, carrier_event[carrier_event$event == 1, ],
                                          sum, drop = FALSE)
    if (nrow(events_by_carrier) < 2L) return("binary_separation")
  }
  if (grepl("Non-finite frailtypack estimate|frailtypack", reason, ignore.case = TRUE)) {
    return("frailtypack_optimizer")
  }
  "unknown"
}
