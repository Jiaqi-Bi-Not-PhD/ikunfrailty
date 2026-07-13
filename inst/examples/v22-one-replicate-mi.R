## V2.2-style one-replicate multiple-imputation example for ikunfrailty.
##
## The data objects loaded below are generated once from the V2.2
## simfam/FamEvent pipeline and shipped with the package. They contain one
## complete selected pop+ replicate with 498 families and three incomplete
## versions made from that same complete data:
##   - ikun_example_continuous_mar20: 20 percent MAR PRS missingness
##   - ikun_example_binary_mar20: 20 percent MAR carrier missingness
##   - ikun_example_joint_mar20: 20 percent MAR PRS and carrier missingness

library(ikunfrailty)
library(survival)

data("ikun_example_continuous_mar20", package = "ikunfrailty")
data("ikun_example_binary_mar20", package = "ikunfrailty")
data("ikun_example_joint_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")
data("ikun_example_pedigree", package = "ikunfrailty")
data("ikun_example_missing_diagnostics", package = "ikunfrailty")
data("ikun_example_metadata", package = "ikunfrailty")

ikun_example_metadata
ikun_example_missing_diagnostics

prior <- pdmi_prior(
  continuous = list(
    newx = normal_kinship(~ mgene, covariance = "kinship+iid")
  ),
  binary = list(
    mgene = carrier_hwe(q = "estimate", q0 = 0.02, n0 = 50)
  )
)

analysis_formula <- Surv(t0, time, status) ~ mgene + newx + cluster(famID)

fit_continuous <- pdmi_frailty(
  formula = analysis_formula,
  data = ikun_example_continuous_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  report = "both",
  pedigree = ikun_example_pedigree,
  seed = 910135,
  progress = TRUE
)

fit_binary <- pdmi_frailty(
  formula = analysis_formula,
  data = ikun_example_binary_mar20,
  kinship = ikun_example_kinship,
  impute = list(binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  report = "both",
  pedigree = ikun_example_pedigree,
  seed = 920135,
  progress = TRUE
)

fit_joint <- pdmi_frailty(
  formula = analysis_formula,
  data = ikun_example_joint_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx", binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  report = "both",
  pedigree = ikun_example_pedigree,
  seed = 930135,
  progress = TRUE
)

summary(fit_continuous)
summary(fit_binary)
summary(fit_joint)

pen_summary(fit_joint)
imputation_model(fit_joint)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- pen_plot(fit_joint)
  print(p)
}
