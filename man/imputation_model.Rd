\name{imputation_model}
\alias{imputation_model}
\alias{imputation_model.pdmi_frailty}
\title{Print the automatically derived congenial imputation model}
\usage{
imputation_model(x, ...)
\method{imputation_model}{pdmi_frailty}(x, ...)
}
\arguments{
\item{x}{A \code{pdmi_frailty} object.}
\item{...}{Unused.}
}
\value{
Invisibly returns a \code{pdmi_imputation_model} object.
}
\description{
Prints the disease likelihood, ascertainment denominator, support factor,
working covariate prior, and sampler pieces derived from the analysis formula
and \code{prior}.
}
