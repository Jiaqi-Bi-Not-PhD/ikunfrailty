test_that("summary, coef, vcov, and penetrance summaries work on a fitted-object shell", {
  omega <- stats::setNames(c(0.8, 4.7, 2.2, 1.0, 0.5), ikunfrailty:::v22_omega_names())
  V <- diag(c(0.01, 0.01, 0.04, 0.04, 0.02))
  dimnames(V) <- list(names(omega), names(omega))
  obj <- list(
    call = quote(pdmi_frailty()),
    report = "both",
    pen_ci = TRUE,
    conf.level = 0.95,
    M = 2L,
    B = 0L,
    numit = 1L,
    missing_type = "continuous",
    prior_version = "C-R",
    result = list(
      convergence = TRUE,
      pooled = list(omega = omega, vcov_omega = V),
      penetrance = data.frame(age = 40, prs = 0, gene = 1,
                              estimate = 0.25, se = 0.05)
    ),
    convergence = TRUE,
    failure_reason = NA_character_,
    diagnostics = list(),
    model = structure(list(
      automatic = TRUE,
      technical_note_target = "test",
      missing_type = "continuous",
      prior_version = "C-R",
      analysis_formula = "Surv(t0, time, status) ~ mgene + newx + cluster(famID)",
      kernel = "kernel",
      disease_component = "disease",
      ascertainment_component = "alpha",
      support_component = "support",
      sampler = list(disease_parameter_draw = "theta",
                     frailty_draw = "frailty",
                     continuous_missing_draw = "continuous",
                     binary_missing_draw = "not active",
                     joint_extra_term = "not active")
    ), class = "pdmi_imputation_model")
  )
  class(obj) <- "pdmi_frailty"

  expect_equal(coef(obj), omega)
  expect_equal(vcov(obj), V)
  expect_s3_class(summary(obj), "summary.pdmi_frailty")
  pen <- pen_summary(obj)
  expect_s3_class(pen, "pdmi_penetrance_summary")
  expect_true(all(c("lower", "upper") %in% names(pen)))
  expect_gte(pen$lower, 0)
  expect_lte(pen$upper, 1)

  skip_if_not_installed("ggplot2")
  expect_s3_class(pen_plot(obj), "ggplot")
})
