#' ikunfrailty: imputation for kinship-induced frailty under nonrandom ascertainment
#'
#' `ikunfrailty` stands for Imputation of covariates for Kinship-induced
#' Unobserved frailty under Nonrandom ascertainment. It implements the proposed
#' posterior-draw multiple imputation (PDMI) engine developed for the V2.2
#' ascertainment-corrected frailtypack technical-note model. The package
#' derives the congenial imputation kernel from the user analysis formula and
#' the working covariate prior:
#'
#' \deqn{
#' f(X_M \mid X_O,Y,P,A=1;\psi) \propto
#' \prod_i
#' \frac{s_i(Y_i,X_i)M_i(Y_i\mid X_i;\theta)}
#' {\alpha_i^F(X_i;\theta)}
#' g_i^A(X_i\mid P_i;\eta_X).
#' }
#'
#' The validated first engine supports the scalar V2.2 setting: one continuous
#' PRS covariate and one binary carrier covariate, separately or jointly
#' missing, in an ascertainment-corrected Weibull/lognormal correlated frailty
#' model fit by `frailtypack::frailtyPenal()`.
#'
#' @seealso [pdmi_frailty()], [pdmi_prior()], [normal_kinship()],
#'   [carrier_hwe()], [pen_summary()], [pen_plot()], [imputation_model()]
#' @aliases ikunfrailty
"_PACKAGE"

utils::globalVariables(c("age", "estimate", "gene", "lower", "profile", "prs", "upper"))
