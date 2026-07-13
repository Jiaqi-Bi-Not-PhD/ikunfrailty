# ikunfrailty

`ikunfrailty` stands for **Imputation of covariates for Kinship-induced
Unobserved frailty under Nonrandom ascertainment**.

This beta research package implements the proposed posterior-draw multiple
imputation method for the V2.2 ascertainment-corrected correlated frailty
model. The current validated engine supports the scalar V2.2 setting:

- one continuous polygenic risk score covariate, such as `newx`;
- one binary carrier covariate, such as `mgene`;
- continuous-only, binary-only, or jointly missing covariates;
- a frailtypack-style ascertainment-corrected correlated frailty analysis
  model with family kinship.

The package derives the congenial imputation kernel from the analysis formula
and the working covariate prior. Users specify the analysis model, data,
kinship matrix, covariates subject to imputation, and prior; `ikunfrailty`
handles the V2.2 posterior-draw imputation engine, completed-data model fits,
Rubin pooling, penetrance summaries, and penetrance plots.

## Installation

Install the beta version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Jiaqi-Bi-Not-PhD/ikunfrailty")
```

To build the vignettes during installation:

```r
remotes::install_github(
  "Jiaqi-Bi-Not-PhD/ikunfrailty",
  build_vignettes = TRUE
)
```

Load the package:

```r
library(ikunfrailty)
```

The package imports `frailtypack`, `survival`, and base R packages. The
penetrance plotting helper uses `ggplot2`; install it if you want plots:

```r
install.packages("ggplot2")
```

## Example Data

The package includes one V2.2-style simulated replicate with 498 selected
families and three incomplete versions:

```r
data("ikun_example_continuous_mar20", package = "ikunfrailty")
data("ikun_example_binary_mar20", package = "ikunfrailty")
data("ikun_example_joint_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")
data("ikun_example_pedigree", package = "ikunfrailty")
```

The incomplete data sets are:

| Data object | Missing covariate pattern |
| --- | --- |
| `ikun_example_continuous_mar20` | `newx` missing, `mgene` observed |
| `ikun_example_binary_mar20` | `mgene` missing, `newx` observed |
| `ikun_example_joint_mar20` | both `newx` and `mgene` missing |

The analysis data include the columns used below:

```r
t0
time
status
mgene
newx
famID
proband
currentage
indID
fatherID
motherID
```

## Specify The Working Prior

For the current scalar V2.2 engine, use `pdmi_prior()` with a continuous PRS
prior and a binary carrier prior:

```r
prior <- pdmi_prior(
  continuous = list(
    newx = normal_kinship(~ mgene, covariance = "kinship+iid")
  ),
  binary = list(
    mgene = carrier_hwe(q = "estimate", q0 = 0.02, n0 = 50)
  )
)
```

The disease part of the imputation model is derived automatically from the
analysis formula. The user supplies only the working covariate prior.

## Quick Smoke Test

The following run uses very small `M`, `B`, and `numit` values only to check
that the installation works. Increase them for real analyses.

```r
library(ikunfrailty)

data("ikun_example_joint_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")
data("ikun_example_pedigree", package = "ikunfrailty")

prior <- pdmi_prior(
  continuous = list(
    newx = normal_kinship(~ mgene, covariance = "kinship+iid")
  ),
  binary = list(
    mgene = carrier_hwe(q = "estimate", q0 = 0.02, n0 = 50)
  )
)

fit <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_joint_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx", binary = "mgene"),
  prior = prior,
  M = 2,
  B = 2,
  numit = 1,
  pedigree = ikun_example_pedigree,
  seed = 930135,
  progress = TRUE
)
```

For a real beta test run, use larger values such as `M = 20`, `B = 50`, and
`numit = 10`, matching the V2.2 examples.

## Parameter Estimates

Use `summary()` for pooled parameter estimates:

```r
summary(fit)
coef(fit)
vcov(fit)
```

Inspect the automatically derived imputation model:

```r
imputation_model(fit)
```

## Penetrance Estimates

Use `pen_summary()` for penetrance estimates and confidence intervals:

```r
pen_summary(fit)
```

Use `pen_plot()` for penetrance curves:

```r
p <- pen_plot(fit)
p
```

To request penetrance only, parameter estimates only, or both:

```r
fit <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_joint_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx", binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  pedigree = ikun_example_pedigree,
  report = "both",
  pen_ci = TRUE,
  progress = TRUE
)
```

## Continuous-Only Missingness

Use `ikun_example_continuous_mar20` when only the continuous PRS covariate
`newx` is missing:

```r
data("ikun_example_continuous_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")

fit_cont <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_continuous_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  progress = TRUE
)

summary(fit_cont)
pen_summary(fit_cont)
```

## Binary-Only Missingness

Use `ikun_example_binary_mar20` when only the binary carrier covariate `mgene`
is missing. Pass the pedigree object for binary or joint imputation:

```r
data("ikun_example_binary_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")
data("ikun_example_pedigree", package = "ikunfrailty")

fit_bin <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_binary_mar20,
  kinship = ikun_example_kinship,
  impute = list(binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  pedigree = ikun_example_pedigree,
  progress = TRUE
)

summary(fit_bin)
pen_summary(fit_bin)
```

## Joint Missingness

Use `ikun_example_joint_mar20` when both `newx` and `mgene` are missing:

```r
data("ikun_example_joint_mar20", package = "ikunfrailty")
data("ikun_example_kinship", package = "ikunfrailty")
data("ikun_example_pedigree", package = "ikunfrailty")

fit_joint <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_joint_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx", binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  pedigree = ikun_example_pedigree,
  progress = TRUE
)

summary(fit_joint)
pen_summary(fit_joint)
pen_plot(fit_joint)
```

If `impute` lists a covariate that has no missing values in `data`, the package
warns and treats that covariate as observed. If a covariate has missing values
but is not listed in `impute`, the package stops with an error.

## Custom Penetrance Grid

By default, penetrance is evaluated over the V2.2-style grid. You can supply a
custom grid:

```r
fit_grid <- pdmi_frailty(
  survival::Surv(t0, time, status) ~ mgene + newx + survival::cluster(famID),
  data = ikun_example_joint_mar20,
  kinship = ikun_example_kinship,
  impute = list(continuous = "newx", binary = "mgene"),
  prior = prior,
  M = 20,
  B = 50,
  numit = 10,
  pedigree = ikun_example_pedigree,
  pen_grid = list(
    age = c(40, 50, 60, 70, 80),
    prs = c(-0.5, 0, 0.5),
    gene = c(0, 1),
    k0 = 1
  ),
  pen_ci = TRUE,
  progress = TRUE
)
```

## Documentation

Package help:

```r
?ikunfrailty
?pdmi_frailty
?pdmi_prior
?pen_summary
?pen_plot
```

Vignettes:

```r
browseVignettes("ikunfrailty")
```

Included vignettes cover:

- a full `pdmi_frailty()` workflow;
- the derived congenial imputation model;
- prior helper functions and arguments;
- one V2.2-style simulated replicate example.

## Beta Scope And Limitations

This beta release is intended for testing the validated V2.2 scalar engine. It
does not yet support arbitrary numbers of continuous or binary covariates,
transformed covariates, interactions, or non-target frailtypack model classes.

The current target model assumes the V2.2 technical-note setting, including:

- Weibull baseline;
- lognormal correlated frailty;
- kinship matrix supplied by the user;
- pop+ ascertainment support;
- observed carrier status for affected carrier probands;
- one continuous PRS covariate and one binary carrier covariate.

Please report installation issues, model failures, or unclear documentation
through GitHub issues.
