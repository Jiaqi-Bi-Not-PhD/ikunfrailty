make_parser_data <- function() {
  data.frame(
    indID = paste0("id", 1:4),
    famID = c(1, 1, 2, 2),
    fatherID = NA_character_,
    motherID = NA_character_,
    t0 = 0,
    time = c(55, 60, 58, 62),
    status = c(1, 0, 1, 0),
    mgene = c(1, 0, 1, 0),
    newx = c(NA, 0.15, 0.20, -0.10),
    proband = c(1, 0, 1, 0),
    currentage = c(55, 60, 58, 62)
  )
}

make_parser_K <- function(dat) {
  K <- diag(nrow(dat))
  rownames(K) <- colnames(K) <- dat$indID
  K
}

test_that("formula parser maps user columns to V2.2 engine columns", {
  dat <- make_parser_data()
  K <- make_parser_K(dat)
  prior <- pdmi_prior(
    continuous = list(newx = normal_kinship(~ mgene, covariance = "kinship+iid"))
  )

  parsed <- ikunfrailty:::pdmi_parse_model_spec(
    survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
    data = dat,
    kinship = K,
    impute = "newx",
    prior = prior,
    id = "indID",
    proband = "proband",
    currentage = "currentage",
    pedigree = NULL
  )

  expect_equal(parsed$spec$missing_type, "continuous")
  expect_equal(parsed$spec$continuous, "newx")
  expect_equal(parsed$spec$binary, "mgene")
  expect_equal(rownames(parsed$K), dat$indID)
  expect_true(all(c("t0", "time", "status", "mgene", "newx", "famID",
                    "proband", "currentage", "indID") %in% names(parsed$dat)))
})

test_that("binary proband missingness is rejected", {
  dat <- make_parser_data()
  dat$newx[1] <- 0.05
  dat$mgene[1] <- NA
  K <- make_parser_K(dat)
  prior <- pdmi_prior(binary = list(mgene = carrier_hwe(q = "estimate")))

  expect_error(
    ikunfrailty:::pdmi_parse_model_spec(
      survival::Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = dat,
      kinship = K,
      impute = "mgene",
      prior = prior,
      id = "indID",
      proband = "proband",
      currentage = "currentage",
      pedigree = NULL
    ),
    "proband carrier"
  )
})

test_that("non-imputed continuous missingness is rejected", {
  dat <- make_parser_data()
  K <- make_parser_K(dat)
  prior <- pdmi_prior(binary = list(mgene = carrier_hwe(q = "estimate")))

  expect_error(
    ikunfrailty:::pdmi_parse_model_spec(
      survival::Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = dat,
      kinship = K,
      impute = "mgene",
      prior = prior,
      id = "indID",
      proband = "proband",
      currentage = "currentage",
      pedigree = NULL
    ),
    "Continuous covariate `newx` contains missing values"
  )
})

test_that("overlisted binary impute target is dropped when fully observed", {
  dat <- make_parser_data()
  K <- make_parser_K(dat)
  prior <- pdmi_prior(
    continuous = list(newx = normal_kinship(~ mgene, covariance = "kinship+iid")),
    binary = list(mgene = carrier_hwe(q = "estimate"))
  )

  parsed <- NULL
  expect_warning(
    parsed <- ikunfrailty:::pdmi_parse_model_spec(
      survival::Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = dat,
      kinship = K,
      impute = list(continuous = "newx", binary = "mgene"),
      prior = prior,
      id = "indID",
      proband = "proband",
      currentage = "currentage",
      pedigree = NULL
    ),
    "Binary covariate `mgene` was listed in `impute` but has no missing values"
  )

  expect_equal(parsed$spec$missing_type, "continuous")
  expect_equal(parsed$spec$impute_vars, "newx")
  expect_equal(parsed$spec$dropped_impute_vars, "mgene")
})

test_that("overlisted continuous impute target is dropped when fully observed", {
  dat <- make_parser_data()
  dat$newx[1] <- 0.05
  dat$mgene[2] <- NA
  K <- make_parser_K(dat)
  prior <- pdmi_prior(
    continuous = list(newx = normal_kinship(~ mgene, covariance = "kinship+iid")),
    binary = list(mgene = carrier_hwe(q = "estimate"))
  )

  parsed <- NULL
  expect_warning(
    parsed <- ikunfrailty:::pdmi_parse_model_spec(
      survival::Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = dat,
      kinship = K,
      impute = list(continuous = "newx", binary = "mgene"),
      prior = prior,
      id = "indID",
      proband = "proband",
      currentage = "currentage",
      pedigree = NULL
    ),
    "Continuous covariate `newx` was listed in `impute` but has no missing values"
  )

  expect_equal(parsed$spec$missing_type, "binary")
  expect_equal(parsed$spec$impute_vars, "mgene")
  expect_equal(parsed$spec$dropped_impute_vars, "newx")
})
