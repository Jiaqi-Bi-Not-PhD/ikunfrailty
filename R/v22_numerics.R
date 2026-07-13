## ============================================================
## Shared numerical utilities.
## ============================================================

v22_safe_exp <- function(x, lim = 700) exp(pmax(pmin(x, lim), -lim))

v22_inv_logit <- function(x) 1 / (1 + exp(-pmax(pmin(x, 35), -35)))

v22_log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}

v22_safe_chol <- function(A, jitter = 1e-8, max_tries = 8L) {
  A <- as.matrix(A)
  A <- 0.5 * (A + t(A))
  for (k in seq_len(max_tries)) {
    R <- tryCatch(chol(A + diag(jitter * 10^(k - 1L), nrow(A))), error = function(e) NULL)
    if (!is.null(R)) return(R)
  }
  ee <- eigen(A, symmetric = TRUE)
  lam <- pmax(ee$values, jitter)
  chol(ee$vectors %*% diag(lam, length(lam)) %*% t(ee$vectors))
}

v22_solve_spd <- function(A, b = NULL, jitter = 1e-8) {
  if (nrow(A) == 1L) {
    val <- max(as.numeric(A[1, 1]), jitter)
    if (is.null(b)) return(matrix(1 / val, 1L, 1L))
    return(as.numeric(b) / val)
  }
  R <- v22_safe_chol(A, jitter = jitter)
  if (is.null(b)) return(chol2inv(R))
  backsolve(R, forwardsolve(t(R), b))
}

v22_rmvnorm_cov <- function(mean, Sigma) {
  mean <- as.numeric(mean)
  if (!length(mean)) return(numeric(0))
  if (length(mean) == 1L) return(mean + rnorm(1L) * sqrt(max(as.numeric(Sigma[1, 1]), 1e-12)))
  R <- v22_safe_chol(Sigma)
  mean + as.numeric(t(R) %*% rnorm(length(mean)))
}

v22_rmvnorm_precision <- function(mean, Q) {
  mean <- as.numeric(mean)
  if (!length(mean)) return(numeric(0))
  if (length(mean) == 1L) return(mean + rnorm(1L) / sqrt(max(as.numeric(Q[1, 1]), 1e-12)))
  R <- v22_safe_chol(Q)
  mean + as.numeric(backsolve(R, rnorm(length(mean))))
}

v22_is_psd <- function(A, tol = 1e-7) {
  if (is.null(A) || any(!is.finite(A))) return(FALSE)
  ev <- tryCatch(eigen(0.5 * (A + t(A)), symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  all(is.finite(ev)) && min(ev) >= -tol
}

v22_near_psd <- function(A, eps = 1e-8) {
  A <- 0.5 * (A + t(A))
  ee <- eigen(A, symmetric = TRUE)
  lam <- pmax(ee$values, eps)
  out <- ee$vectors %*% diag(lam, length(lam)) %*% t(ee$vectors)
  dimnames(out) <- dimnames(A)
  0.5 * (out + t(out))
}

v22_matrix_sqrt_sym <- function(A, inverse = FALSE, eps = 1e-10) {
  A <- 0.5 * (A + t(A))
  ee <- eigen(A, symmetric = TRUE)
  lam <- pmax(ee$values, eps)
  lam <- if (isTRUE(inverse)) 1 / sqrt(lam) else sqrt(lam)
  ee$vectors %*% diag(lam, length(lam)) %*% t(ee$vectors)
}

v22_gh_rule <- function(order = 20L) {
  if (order == 20L) {
    nodes <- c(-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,
               -2.78880605843,-2.25497400209,-1.73853771212,-1.23407621540,
               -0.737473728545,-0.245340708301, 0.245340708301, 0.737473728545,
               1.23407621540, 1.73853771212, 2.25497400209, 2.78880605843,
               3.34785456738, 3.94476404012, 4.60368244955, 5.38748089001)
    weights <- c(2.22939364553415129252e-13, 4.39934099227318055360e-10,
                 1.08606937076928169400e-07, 7.80255647853206369415e-06,
                 2.28338636016353967257e-04, 3.24377334223786183218e-03,
                 2.48105208874636108822e-02, 1.09017206020023320014e-01,
                 2.86675505362834129720e-01, 4.62243669600610089650e-01,
                 4.62243669600610089650e-01, 2.86675505362834129720e-01,
                 1.09017206020023320014e-01, 2.48105208874636108822e-02,
                 3.24377334223786183218e-03, 2.28338636016353967257e-04,
                 7.80255647853206369410e-06, 1.08606937076928169400e-07,
                 4.39934099227318055363e-10, 2.22939364553415129252e-13)
    return(list(nodes = nodes, weights = weights))
  }
  if (order == 10L) {
    nodes <- c(-3.43615911884,-2.53273167423,-1.75668364930,-1.03661082979,
               -0.342901327224, 0.342901327224, 1.03661082979, 1.75668364930,
               2.53273167423, 3.43615911884)
    weights <- c(7.64043285523e-06, 1.34364574678e-03, 3.38743944555e-02,
                 2.40138611082e-01, 6.10862633735e-01, 6.10862633735e-01,
                 2.40138611082e-01, 3.38743944555e-02, 1.34364574678e-03,
                 7.64043285523e-06)
    return(list(nodes = nodes, weights = weights))
  }
  stop("v22_gh_rule() currently supports order 10 or 20.")
}

v22_H0 <- function(t, omega, agemin = 0) {
  th <- v22_theta_from_omega(omega)
  td <- pmax(as.numeric(t) - agemin, 0)
  (td / th$lambda)^th$rho
}

v22_H0_diff <- function(time, t0, omega, agemin = 0) {
  pmax(0, v22_H0(time, omega, agemin) - v22_H0(t0, omega, agemin))
}

v22_log_h0 <- function(time, omega, agemin = 0) {
  th <- v22_theta_from_omega(omega)
  td <- pmax(as.numeric(time) - agemin, 1e-12)
  log(th$rho) - th$log_lambda + (th$rho - 1) * (log(td) - th$log_lambda)
}

v22_alpha_popplus <- function(age, x_c, x_b, k_diag, omega, agemin = 0, gh_order = 20L,
                             eps = 1e-12) {
  th <- v22_theta_from_omega(omega)
  gh <- v22_gh_rule(gh_order)
  Hc <- v22_H0(age, omega, agemin)
  v <- pmax(th$sigma_u2 * k_diag, 1e-12)
  eta <- th$beta_c * x_c + th$beta_b * x_b
  d <- Hc * v22_safe_exp(eta + sqrt(2 * v) * gh$nodes)
  S <- sum((gh$weights / sqrt(pi)) * exp(-pmin(d, 745)))
  pmax(pmin(1 - S, 1 - eps), eps)
}

v22_penetrance <- function(age, x_c, x_b, omega, k0 = 1, agemin = 0, gh_order = 20L) {
  th <- v22_theta_from_omega(omega)
  gh <- v22_gh_rule(gh_order)
  Ht <- v22_H0(age, omega, agemin)
  v <- pmax(th$sigma_u2 * k0, 1e-12)
  eta <- th$beta_c * x_c + th$beta_b * x_b
  d <- Ht * v22_safe_exp(eta + sqrt(2 * v) * gh$nodes)
  1 - sum((gh$weights / sqrt(pi)) * exp(-pmin(d, 745)))
}

v22_penetrance_grid <- function(omega, config = v22_default_config()) {
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
  grid
}

v22_penetrance_gradient_omega <- function(age, x_c, x_b, omega, k0 = 1, agemin = 0,
                                         gh_order = 20L) {
  th <- v22_theta_from_omega(omega)
  gh <- v22_gh_rule(gh_order)
  td <- pmax(age - agemin, 0)
  if (td <= 0) return(setNames(rep(0, 5), v22_omega_names()))
  Ht <- v22_H0(age, omega, agemin)
  v <- pmax(th$sigma_u2 * k0, 1e-12)
  q <- sqrt(2 * v) * gh$nodes
  d <- Ht * v22_safe_exp(th$beta_c * x_c + th$beta_b * x_b + q)
  common <- (gh$weights / sqrt(pi)) * exp(-pmin(d, 745)) * d
  qprime_sigma <- gh$nodes * k0 / sqrt(2 * v)
  grad_theta <- c(
    beta_c = sum(common * x_c),
    beta_b = sum(common * x_b),
    lambda = sum(common * (-th$rho / th$lambda)),
    rho = sum(common * log(td / th$lambda)),
    sigma_u2 = sum(common * qprime_sigma)
  )
  out <- c(
    log.rho = grad_theta["rho"] * th$rho,
    log.lambda = grad_theta["lambda"] * th$lambda,
    beta_b = grad_theta["beta_b"],
    beta_c = grad_theta["beta_c"],
    sigma_u2 = grad_theta["sigma_u2"]
  )
  setNames(as.numeric(out), v22_omega_names())
}

v22_slice_univariate <- function(f, x0, w = 0.7, m = 80L, lower = -Inf, upper = Inf) {
  fx0 <- f(x0)
  if (!is.finite(fx0)) {
    x0 <- 0
    fx0 <- f(x0)
  }
  if (!is.finite(fx0)) stop("Slice sampler was initialized at a non-finite log density.")
  y <- fx0 - rexp(1L)
  L <- x0 - runif(1L, 0, w)
  R <- L + w
  J <- floor(runif(1L, 0, m))
  K <- (m - 1L) - J
  while (J > 0 && L > lower && is.finite(f(L)) && f(L) > y) {
    L <- L - w
    J <- J - 1L
  }
  while (K > 0 && R < upper && is.finite(f(R)) && f(R) > y) {
    R <- R + w
    K <- K - 1L
  }
  repeat {
    lo <- max(L, lower)
    hi <- min(R, upper)
    if (!(lo < hi)) return(x0)
    x1 <- runif(1L, lo, hi)
    fx1 <- f(x1)
    if (is.finite(fx1) && fx1 >= y) return(x1)
    if (x1 < x0) L <- x1 else R <- x1
  }
}

v22_family_blocks <- function(dat) split(seq_len(nrow(dat)), as.factor(dat$famID))

v22_align_K <- function(K, dat, id_col = "indID") {
  ids <- as.character(dat[[id_col]])
  if (!is.null(rownames(K))) {
    ord <- match(ids, rownames(K))
    if (anyNA(ord)) stop("K rownames do not match dat$", id_col, ".")
    K <- K[ord, ord, drop = FALSE]
  }
  K <- as.matrix(K)
  rownames(K) <- colnames(K) <- ids
  0.5 * (K + t(K))
}

v22_precompute_K_blocks <- function(K, dat, jitter = 1e-8) {
  blocks <- v22_family_blocks(dat)
  lapply(blocks, function(idx) {
    Ki_raw <- as.matrix(K[idx, idx, drop = FALSE])
    Ki <- 0.5 * (Ki_raw + t(Ki_raw))
    R <- v22_safe_chol(Ki, jitter = jitter)
    list(idx = idx, K = Ki, K_inv = chol2inv(R),
         logdetK = 2 * sum(log(diag(R))), n = length(idx))
  })
}

v22_make_kinship_cache <- function(K, dat, id_col = "indID") {
  K <- v22_align_K(K, dat, id_col = id_col)
  ids <- as.character(dat[[id_col]])
  list(
    ids = ids,
    K = K,
    blocks = v22_precompute_K_blocks(K, dat),
    n = length(ids),
    n_families = length(unique(dat$famID))
  )
}

v22_kinship_cache_matches <- function(cache, dat, id_col = "indID") {
  !is.null(cache) &&
    !is.null(cache$ids) &&
    identical(cache$ids, as.character(dat[[id_col]]))
}

v22_get_cached_K <- function(config, K, dat) {
  cache <- config$kinship_cache %||% NULL
  if (v22_kinship_cache_matches(cache, dat)) return(cache$K)
  v22_align_K(K, dat)
}

v22_get_cached_K_blocks <- function(config, K, dat) {
  cache <- config$kinship_cache %||% NULL
  if (v22_kinship_cache_matches(cache, dat)) return(cache$blocks)
  v22_precompute_K_blocks(K, dat)
}

v22_proband_indices <- function(dat) {
  blocks <- v22_family_blocks(dat)
  vapply(blocks, function(idx) {
    loc <- which(as.numeric(dat$proband[idx]) == 1)
    if (length(loc) != 1L) return(NA_integer_)
    idx[loc]
  }, integer(1))
}
