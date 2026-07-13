pdmi_call_name <- function(x) {
  if (!is.call(x)) return("")
  nm <- as.character(x[[1]])
  nm[length(nm)]
}

pdmi_symbol_name <- function(x, label) {
  if (!is.symbol(x)) {
    stop("The ", label, " must be a column name, not an expression.",
         call. = FALSE)
  }
  as.character(x)
}

pdmi_surv_columns <- function(formula) {
  lhs <- formula[[2]]
  if (!is.call(lhs) || !identical(pdmi_call_name(lhs), "Surv")) {
    stop("`formula` must have a survival response such as ",
         "Surv(t0, time, status) or Surv(time, status).", call. = FALSE)
  }
  args <- as.list(lhs)[-1]
  if (length(args) == 2L) {
    list(
      start = NULL,
      time = pdmi_symbol_name(args[[1]], "Surv time argument"),
      status = pdmi_symbol_name(args[[2]], "Surv status argument")
    )
  } else if (length(args) == 3L) {
    list(
      start = pdmi_symbol_name(args[[1]], "Surv start argument"),
      time = pdmi_symbol_name(args[[2]], "Surv stop/time argument"),
      status = pdmi_symbol_name(args[[3]], "Surv status argument")
    )
  } else {
    stop("The V2.2 engine supports Surv(time, status) or ",
         "Surv(t0, time, status).", call. = FALSE)
  }
}

pdmi_formula_terms <- function(formula) {
  trm <- stats::terms(formula, keep.order = TRUE)
  labels <- attr(trm, "term.labels")
  cluster_idx <- grep("(^|::)cluster\\(", labels)
  if (length(cluster_idx) != 1L) {
    stop("`formula` must contain exactly one cluster(family) term.",
         call. = FALSE)
  }
  cluster_term <- labels[cluster_idx]
  cluster_col <- sub("^.*cluster\\((.*)\\)$", "\\1", cluster_term)
  cluster_col <- trimws(cluster_col)
  if (!nzchar(cluster_col) || grepl("[+*:()/]", cluster_col)) {
    stop("The cluster() argument must be a single family column name.",
         call. = FALSE)
  }

  covariates <- labels[-cluster_idx]
  if (length(covariates) != 2L) {
    stop("The validated V2.2 engine requires exactly two disease covariates: ",
         "one continuous PRS and one binary carrier covariate.", call. = FALSE)
  }
  if (any(grepl("[:()*^/+-]", covariates))) {
    stop("The validated V2.2 engine does not support transformed, interaction, ",
         "or arithmetic covariate terms.", call. = FALSE)
  }
  list(covariates = covariates, cluster = cluster_col)
}

pdmi_parse_impute <- function(impute) {
  if (missing(impute) || is.null(impute)) {
    stop("`impute` must name the covariate(s) subject to missingness.",
         call. = FALSE)
  }
  hints <- list(continuous = character(0), binary = character(0))
  if (is.character(impute)) {
    vars <- impute
  } else if (is.list(impute)) {
    vars <- unique(as.character(unlist(impute, use.names = FALSE)))
    nms <- names(impute) %||% rep("", length(impute))
    for (j in seq_along(impute)) {
      nm <- nms[[j]]
      if (nm %in% names(hints)) {
        hints[[nm]] <- unique(c(hints[[nm]], as.character(unlist(impute[[j]]))))
      }
    }
  } else {
    stop("`impute` must be a character vector or a list.", call. = FALSE)
  }
  vars <- unique(vars[nzchar(vars)])
  if (!length(vars)) stop("`impute` did not contain any variable names.",
                          call. = FALSE)
  list(vars = vars, hints = hints)
}

pdmi_is_binary_like <- function(x) {
  if (is.logical(x)) return(TRUE)
  if (is.factor(x)) {
    vals <- stats::na.omit(as.character(x))
    return(length(vals) > 0L && all(vals %in% c("0", "1")))
  }
  vals <- stats::na.omit(as.numeric(x))
  length(vals) > 0L && all(vals %in% c(0, 1))
}

pdmi_as_binary <- function(x, label) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) {
    vals <- as.character(x)
    if (!all(stats::na.omit(vals) %in% c("0", "1"))) {
      stop("Binary covariate `", label, "` must be coded 0/1.", call. = FALSE)
    }
    return(as.numeric(vals))
  }
  vals <- as.numeric(x)
  ok <- is.na(vals) | vals %in% c(0, 1)
  if (!all(ok)) stop("Binary covariate `", label, "` must be coded 0/1.",
                     call. = FALSE)
  vals
}

pdmi_as_numeric_column <- function(x, label) {
  vals <- suppressWarnings(as.numeric(x))
  if (all(is.na(vals)) && any(!is.na(x))) {
    stop("Column `", label, "` must be numeric.", call. = FALSE)
  }
  vals
}

pdmi_resolve_covariates <- function(covariates, data, prior, impute_spec) {
  if (!inherits(prior, "pdmi_prior")) {
    stop("`prior` must be created by pdmi_prior().", call. = FALSE)
  }
  prior_vars <- c(names(prior$continuous), names(prior$binary))
  bad_prior <- setdiff(prior_vars, covariates)
  if (length(bad_prior)) {
    stop("Prior covariate(s) not in the analysis formula: ",
         paste(bad_prior, collapse = ", "), call. = FALSE)
  }

  continuous <- intersect(covariates, names(prior$continuous))
  binary <- intersect(covariates, names(prior$binary))
  continuous <- unique(c(continuous, intersect(covariates, impute_spec$hints$continuous)))
  binary <- unique(c(binary, intersect(covariates, impute_spec$hints$binary)))

  if (length(continuous) > 1L || length(binary) > 1L) {
    stop("The scalar V2.2 engine supports exactly one continuous and one binary covariate.",
         call. = FALSE)
  }

  if (!length(binary)) {
    binary <- covariates[vapply(covariates, function(v) pdmi_is_binary_like(data[[v]]),
                                logical(1))]
  }
  if (length(binary) != 1L) {
    stop("Could not identify the single binary carrier covariate. ",
         "Name it in pdmi_prior(binary = ...) or impute = list(binary = ...).",
         call. = FALSE)
  }

  if (!length(continuous)) continuous <- setdiff(covariates, binary)
  if (length(continuous) != 1L || identical(continuous, binary)) {
    stop("Could not identify the single continuous PRS covariate.",
         call. = FALSE)
  }
  if (pdmi_is_binary_like(data[[continuous]])) {
    stop("Continuous covariate `", continuous, "` appears to be binary.",
         call. = FALSE)
  }

  impute_vars <- impute_spec$vars
  bad_impute <- setdiff(impute_vars, covariates)
  if (length(bad_impute)) {
    stop("`impute` variable(s) not in the analysis formula: ",
         paste(bad_impute, collapse = ", "), call. = FALSE)
  }
  requested_cont <- continuous %in% impute_vars
  requested_bin <- binary %in% impute_vars
  actual_cont_missing <- anyNA(data[[continuous]])
  actual_bin_missing <- anyNA(data[[binary]])

  if (actual_cont_missing && !requested_cont) {
    stop("Continuous covariate `", continuous,
         "` contains missing values but was not listed in `impute`.",
         call. = FALSE)
  }
  if (actual_bin_missing && !requested_bin) {
    stop("Binary covariate `", binary,
         "` contains missing values but was not listed in `impute`.",
         call. = FALSE)
  }
  if (requested_cont && !actual_cont_missing) {
    warning("Continuous covariate `", continuous,
            "` was listed in `impute` but has no missing values; ",
            "it will be treated as observed.", call. = FALSE)
  }
  if (requested_bin && !actual_bin_missing) {
    warning("Binary covariate `", binary,
            "` was listed in `impute` but has no missing values; ",
            "it will be treated as observed.", call. = FALSE)
  }

  cont_missing <- requested_cont && actual_cont_missing
  bin_missing <- requested_bin && actual_bin_missing
  missing_type <- if (cont_missing && bin_missing) {
    "joint"
  } else if (cont_missing) {
    "continuous"
  } else if (bin_missing) {
    "binary"
  } else {
    stop("`impute` must name at least one analysis covariate that contains missing values.",
         call. = FALSE)
  }

  effective_impute_vars <- c(if (cont_missing) continuous else character(0),
                             if (bin_missing) binary else character(0))
  dropped_impute_vars <- setdiff(impute_vars, effective_impute_vars)

  list(continuous = continuous,
       binary = binary,
       missing_type = missing_type,
       impute_vars = effective_impute_vars,
       requested_impute_vars = impute_vars,
       dropped_impute_vars = dropped_impute_vars)
}

pdmi_require_columns <- function(data, cols, context = "`data`") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop(context, " is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

pdmi_optional_column <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) data[[hit[[1]]]] else rep(NA_character_, nrow(data))
}

pdmi_prepare_analysis_data <- function(data, surv, terms, roles,
                                       id, proband, currentage) {
  data <- as.data.frame(data)
  required <- c(surv$time, surv$status, terms$cluster,
                roles$continuous, roles$binary, id, proband, currentage)
  if (!is.null(surv$start)) required <- c(required, surv$start)
  pdmi_require_columns(data, unique(required))

  out <- data.frame(
    t0 = if (is.null(surv$start)) rep(0, nrow(data)) else
      pdmi_as_numeric_column(data[[surv$start]], surv$start),
    time = pdmi_as_numeric_column(data[[surv$time]], surv$time),
    status = pdmi_as_numeric_column(data[[surv$status]], surv$status),
    mgene = pdmi_as_binary(data[[roles$binary]], roles$binary),
    newx = pdmi_as_numeric_column(data[[roles$continuous]], roles$continuous),
    famID = data[[terms$cluster]],
    proband = pdmi_as_numeric_column(data[[proband]], proband),
    currentage = pdmi_as_numeric_column(data[[currentage]], currentage),
    indID = data[[id]],
    stringsAsFactors = FALSE
  )
  out$fatherID <- pdmi_optional_column(data, c("fatherID", "father", "fatherid"))
  out$motherID <- pdmi_optional_column(data, c("motherID", "mother", "motherid"))

  if (any(is.na(out$indID)) || any(!nzchar(as.character(out$indID)))) {
    stop("Individual id column `", id, "` contains missing or blank values.",
         call. = FALSE)
  }
  if (anyDuplicated(as.character(out$indID))) {
    stop("Individual id column `", id, "` must be unique in the analysis data.",
         call. = FALSE)
  }
  if (any(is.na(out$famID))) stop("Family cluster column contains missing values.",
                                  call. = FALSE)
  if (any(!is.finite(out$t0)) || any(!is.finite(out$time)) ||
      any(!is.finite(out$status)) || any(!is.finite(out$proband)) ||
      any(!is.finite(out$currentage))) {
    stop("Survival time/status, proband, and current age columns must be finite.",
         call. = FALSE)
  }
  if (!all(out$status %in% c(0, 1))) {
    stop("The event/status column must be coded 0/1.", call. = FALSE)
  }
  out
}

pdmi_validate_missingness <- function(dat, roles) {
  if (!roles$continuous %in% roles$impute_vars && anyNA(dat$newx)) {
    stop("Continuous covariate `", roles$continuous,
         "` contains missing values but was not listed in `impute`.",
         call. = FALSE)
  }
  if (!roles$binary %in% roles$impute_vars && anyNA(dat$mgene)) {
    stop("Binary covariate `", roles$binary,
         "` contains missing values but was not listed in `impute`.",
         call. = FALSE)
  }
  if (any(is.na(dat$mgene) & as.numeric(dat$proband) == 1)) {
    stop("The current V2.2 binary/joint engine requires proband carrier status ",
         "to be observed and equal to 1 under pop+ ascertainment.", call. = FALSE)
  }
  if (!v22_check_popplus_support(dat)) {
    stop("Analysis data violate pop+ support: each family must have exactly one ",
         "affected carrier proband.", call. = FALSE)
  }
  invisible(TRUE)
}

pdmi_prepare_kinship <- function(kinship, dat) {
  K <- as.matrix(kinship)
  if (!is.numeric(K) || nrow(K) != ncol(K)) {
    stop("`kinship` must be a square numeric matrix.", call. = FALSE)
  }
  ids <- as.character(dat$indID)
  if (nrow(K) != length(ids)) {
    stop("`kinship` dimensions must match the number of analysis rows.",
         call. = FALSE)
  }
  if (is.null(rownames(K))) rownames(K) <- ids
  if (is.null(colnames(K))) colnames(K) <- rownames(K)
  if (!all(ids %in% rownames(K)) || !all(ids %in% colnames(K))) {
    stop("`kinship` row and column names must contain all analysis ids.",
         call. = FALSE)
  }
  K <- K[ids, ids, drop = FALSE]
  0.5 * (K + t(K))
}

pdmi_pick_pedigree_column <- function(pedigree, candidates, default = NULL) {
  hit <- candidates[candidates %in% names(pedigree)]
  if (length(hit)) pedigree[[hit[[1]]]] else default
}

pdmi_prepare_pedigree <- function(pedigree, data, dat, terms, roles, id, proband) {
  if (is.null(pedigree)) return(dat)
  pedigree <- as.data.frame(pedigree)
  ped_id <- pdmi_pick_pedigree_column(pedigree, c(id, "indID", "id"))
  ped_fam <- pdmi_pick_pedigree_column(pedigree, c(terms$cluster, "famID", "family"))
  if (is.null(ped_id) || is.null(ped_fam)) {
    stop("`pedigree` must contain individual id and family id columns.",
         call. = FALSE)
  }
  out <- data.frame(
    indID = ped_id,
    famID = ped_fam,
    fatherID = pdmi_pick_pedigree_column(pedigree, c("fatherID", "father", "fatherid"),
                                         rep(NA_character_, nrow(pedigree))),
    motherID = pdmi_pick_pedigree_column(pedigree, c("motherID", "mother", "motherid"),
                                         rep(NA_character_, nrow(pedigree))),
    proband = pdmi_pick_pedigree_column(pedigree, c(proband, "proband"),
                                        rep(0, nrow(pedigree))),
    mgene = pdmi_pick_pedigree_column(pedigree, c(roles$binary, "mgene"),
                                      rep(NA_real_, nrow(pedigree))),
    stringsAsFactors = FALSE
  )
  match_idx <- match(as.character(dat$indID), as.character(out$indID))
  if (anyNA(match_idx)) {
    stop("Every analysis id must appear in `pedigree`.", call. = FALSE)
  }
  out$proband[match_idx] <- dat$proband
  out$mgene[match_idx] <- dat$mgene
  out
}

pdmi_parse_model_spec <- function(formula, data, kinship, impute, prior,
                                  id, proband, currentage, pedigree) {
  formula <- stats::as.formula(formula)
  surv <- pdmi_surv_columns(formula)
  terms <- pdmi_formula_terms(formula)
  impute_spec <- pdmi_parse_impute(impute)
  roles <- pdmi_resolve_covariates(terms$covariates, data, prior, impute_spec)
  dat <- pdmi_prepare_analysis_data(data, surv, terms, roles, id, proband, currentage)
  pdmi_validate_missingness(dat, roles)
  K <- pdmi_prepare_kinship(kinship, dat)
  ped <- pdmi_prepare_pedigree(pedigree, data, dat, terms, roles, id, proband)

  spec <- c(
    list(
      formula = formula,
      surv = surv,
      cluster = terms$cluster,
      id = id,
      proband = proband,
      currentage = currentage,
      analysis_columns = list(
        t0 = surv$start %||% "0",
        time = surv$time,
        status = surv$status,
        family = terms$cluster,
        continuous = roles$continuous,
        binary = roles$binary,
        id = id,
        proband = proband,
        currentage = currentage
      ),
      engine_columns = list(
        t0 = "t0",
        time = "time",
        status = "status",
        family = "famID",
        continuous = "newx",
        binary = "mgene",
        id = "indID",
        proband = "proband",
        currentage = "currentage"
      )
    ),
    roles
  )
  list(spec = spec, dat = dat, K = K, pedigree = ped)
}

pdmi_build_kernel_model <- function(spec, prior, prior_version, B, numit) {
  continuous_prior <- prior$continuous[[spec$continuous]]
  binary_prior <- prior$binary[[spec$binary]]
  out <- list(
    automatic = TRUE,
    technical_note_target = "V2.2 scalar ascertainment-corrected frailtypack PDMI",
    missing_type = spec$missing_type,
    prior_version = prior_version,
    analysis_formula = paste(deparse(spec$formula), collapse = " "),
    frailtypack_model = paste(
      "Surv(t0, time, status) ~ mgene + newx + cluster(famID);",
      "hazard = Weibull; RandDist = LogN; recurrentAG = TRUE;",
      "covMatrix1 = kinship; proband/currentage ascertainment correction"
    ),
    kernel = paste(
      "prod_i s_i(Y_i,X_i) M_i(Y_i | X_i; theta) /",
      "alpha_i^F(X_i; theta) * g_i^A(X_i | P_i; eta_X)"
    ),
    disease_component = paste(
      "M_i is the Weibull baseline, lognormal correlated-frailty family",
      "likelihood induced by the analysis formula."
    ),
    ascertainment_component = paste(
      "alpha_i^F is the pop+ proband ascertainment probability computed from",
      "currentage, proband carrier support, PRS, carrier status, kinship diagonal,",
      "and theta by Gauss-Hermite quadrature."
    ),
    support_component = paste(
      "s_i enforces one affected carrier proband per family; binary proband",
      "carrier status is fixed observed at 1."
    ),
    working_covariate_prior = list(
      continuous = continuous_prior,
      binary = binary_prior
    ),
    sampler = list(
      disease_parameter_draw = "exact univariate slice on omega with frailty-marginal Laplace target",
      frailty_draw = "elliptical slice within family",
      continuous_missing_draw = if (spec$missing_type %in% c("continuous", "joint"))
        "elliptical slice from conditional normal prior times disease and -log alpha when proband PRS is missing" else "not active",
      binary_missing_draw = if (spec$missing_type %in% c("binary", "joint"))
        "Mendelian genotype Gibbs update with carrier projection; alpha constant for binary nonproband updates" else "not active",
      joint_extra_term = if (spec$missing_type == "joint")
        "binary genotype update includes the continuous prior density when the PRS prior depends on carrier status" else "not active",
      burn_in = B,
      retained_update_iterations = numit
    )
  )
  class(out) <- "pdmi_imputation_model"
  out
}
