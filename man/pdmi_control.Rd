\name{pdmi_control}
\alias{pdmi_control}
\title{Control settings for pdmi_frailty}
\description{
Build a control list for numerical and sampler options used by
\code{\link{pdmi_frailty}}.
}
\usage{
pdmi_control(
  agemin = 0,
  frailtypack_maxit = 35L,
  gh_order = 20L,
  theta_slice_widths = NULL,
  theta_slice_sweeps = 1L,
  theta_slice_m = 80L,
  continuous_slice_step = 0.7,
  continuous_slice_m = 80L,
  binary_genotype_sweeps = 8L,
  binary_prior_gibbs_burnin = 30L,
  binary_prior_gibbs_draws = 20L,
  binary_prior_gibbs_thin = 2L,
  pdmi_debug = FALSE,
  ncores = 1L,
  ...
)
}
\arguments{
\item{agemin}{Minimum analysis age.}
\item{frailtypack_maxit}{Maximum frailtypack iterations.}
\item{gh_order}{Gauss-Hermite quadrature order, 10 or 20.}
\item{theta_slice_widths}{Named slice widths for disease parameters.}
\item{theta_slice_sweeps}{Number of disease-parameter slice sweeps.}
\item{theta_slice_m}{Maximum stepping-out steps for disease parameters.}
\item{continuous_slice_step}{Continuous sampler step size retained for
compatibility.}
\item{continuous_slice_m}{Maximum continuous-sampler stepping-out steps.}
\item{binary_genotype_sweeps}{Number of genotype Gibbs sweeps per update.}
\item{binary_prior_gibbs_burnin}{Burn-in for binary prior routines.}
\item{binary_prior_gibbs_draws}{Retained binary prior Gibbs draws.}
\item{binary_prior_gibbs_thin}{Binary prior thinning interval.}
\item{pdmi_debug}{Whether to keep short diagnostic traces.}
\item{ncores}{Reserved for future parallel completed-data fits.}
\item{...}{Additional V2.2 configuration overrides.}
}
\value{A \code{pdmi_control} list.}
