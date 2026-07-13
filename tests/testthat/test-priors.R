test_that("prior constructors classify supported V2.2 prior versions", {
  prior <- pdmi_prior(
    continuous = list(newx = normal_kinship(~ mgene, covariance = "kinship+iid")),
    binary = list(mgene = carrier_hwe(q = "estimate"))
  )
  expect_s3_class(prior, "pdmi_prior")

  expect_equal(ikunfrailty:::pdmi_prior_version_continuous(prior$continuous$newx), "C-R")
  expect_equal(ikunfrailty:::pdmi_prior_version_binary(prior$binary$mgene), "B-R")

  oracle <- pdmi_prior(
    continuous = list(newx = normal_kinship(~ 1, covariance = "kinship",
                                            sigma2 = 0.1, estimate = FALSE)),
    binary = list(mgene = carrier_hwe(q = 0.02))
  )
  expect_equal(ikunfrailty:::pdmi_prior_version_continuous(oracle$continuous$newx), "C-O")
  expect_equal(ikunfrailty:::pdmi_prior_version_binary(oracle$binary$mgene), "B-O")
})

test_that("unsupported prior placeholders are not accepted by the scalar engine", {
  bad_cont <- mvn_kinship(~ mgene)
  expect_error(ikunfrailty:::pdmi_prior_version_continuous(bad_cont),
               "normal_kinship")

  bad_bin <- bernoulli_glm(mgene ~ newx)
  expect_error(ikunfrailty:::pdmi_prior_version_binary(bad_bin),
               "carrier_hwe")
})
