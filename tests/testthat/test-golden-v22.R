test_that("public wrapper matches fixed V2.2 golden outputs", {
  or_else <- function(a, b) if (!is.null(a)) a else b
  skip_if_not(
    identical(Sys.getenv("PDMI_RUN_GOLDEN_TESTS"), "true"),
    "Set PDMI_RUN_GOLDEN_TESTS=true and provide V2.2 fixture outputs to run golden tests."
  )
  fixture_path <- testthat::test_path("fixtures", "v22_scalar_golden.rds")
  skip_if_not(file.exists(fixture_path), "Missing tests/testthat/fixtures/v22_scalar_golden.rds")
  skip_if_not_installed("frailtypack")

  fixtures <- readRDS(fixture_path)
  for (scenario in c("continuous", "binary", "joint")) {
    fx <- fixtures[[scenario]]
    fit <- pdmi_frailty(
      formula = survival::Surv(t0, time, status) ~ mgene + newx + cluster(famID),
      data = fx$data,
      kinship = fx$kinship,
      impute = fx$impute,
      prior = fx$prior,
      M = fx$M,
      B = fx$B,
      numit = fx$numit,
      pedigree = fx$pedigree,
      seed = fx$seed,
      control = or_else(fx$control, pdmi_control())
    )
    expect_equal(coef(fit), fx$coef, tolerance = or_else(fx$tolerance, 1e-6))
    expect_equal(pen_summary(fit)$estimate, fx$penetrance$estimate,
                 tolerance = or_else(fx$tolerance, 1e-6))
  }
})
