test_that("packaged V2.2 example data are aligned and incomplete as documented", {
  data("ikun_example_complete", package = "ikunfrailty")
  data("ikun_example_continuous_mar20", package = "ikunfrailty")
  data("ikun_example_binary_mar20", package = "ikunfrailty")
  data("ikun_example_joint_mar20", package = "ikunfrailty")
  data("ikun_example_kinship", package = "ikunfrailty")
  data("ikun_example_pedigree", package = "ikunfrailty")
  data("ikun_example_missing_diagnostics", package = "ikunfrailty")

  expect_equal(length(unique(ikun_example_complete$famID)), 498L)
  expect_equal(nrow(ikun_example_continuous_mar20), nrow(ikun_example_complete))
  expect_equal(nrow(ikun_example_binary_mar20), nrow(ikun_example_complete))
  expect_equal(nrow(ikun_example_joint_mar20), nrow(ikun_example_complete))
  expect_equal(nrow(ikun_example_kinship), nrow(ikun_example_complete))
  expect_equal(rownames(ikun_example_kinship), as.character(ikun_example_complete$indID))
  expect_equal(colnames(ikun_example_kinship), as.character(ikun_example_complete$indID))
  expect_true(all(ikun_example_complete$indID %in% ikun_example_pedigree$indID))

  expect_equal(sum(is.na(ikun_example_continuous_mar20$newx)), 585L)
  expect_equal(sum(is.na(ikun_example_binary_mar20$mgene)), 485L)
  expect_equal(sum(is.na(ikun_example_joint_mar20$newx)), 585L)
  expect_equal(sum(is.na(ikun_example_joint_mar20$mgene)), 485L)

  expect_true(all(!is.na(ikun_example_binary_mar20$mgene[
    ikun_example_binary_mar20$proband == 1
  ])))
  expect_true(all(!is.na(ikun_example_joint_mar20$mgene[
    ikun_example_joint_mar20$proband == 1
  ])))
  expect_equal(ikun_example_missing_diagnostics$continuous$realized_missing_rate, 0.2)
  expect_equal(ikun_example_missing_diagnostics$joint$realized_missing_rate_x, 0.2)

  cache <- ikunfrailty:::v22_make_kinship_cache(ikun_example_kinship, ikun_example_complete)
  expect_true(ikunfrailty:::v22_kinship_cache_matches(cache, ikun_example_complete))
  expect_equal(length(cache$blocks), 498L)
})
