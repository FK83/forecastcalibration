% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_procs.R
\name{knueppel_test}
\alias{knueppel_test}
\title{Knueppel (2015) test for calibration of multi-step forecast distribution}
\usage{
knueppel_test(x, lags = NULL, prewhite = FALSE, meth = "qs", moments = "1234")
}
\arguments{
\item{x}{series of PIT values. Under the null of calibration, unconditional
distribution of \code{x} is standard uniform.}

\item{lags, prewhite, meth}{Parameters for computing the autocorrelation-robust variance
that enters the test statistic.}

\item{moments}{Choice of moments of uniform distribution to be considered. String, either \code{"1234"} (first four moments, the default) or \code{"12"} (first two moments).}
}
\description{
Knueppel (2015) test for calibration of multi-step forecast distribution
}
\details{
The variants of the test are called "$alpha_{1234}^{0}$ and "$alpha_{12}^{0}$ in Knueppel (2015).
}
\examples{
set.seed(1)

# Simulate data that satisfies the null
u1 <- pnorm(arima.sim(model = list(ar = .5), n = 2000, sd = sqrt(1-.5^2)))
# Show PIT histogram
hist(u1)
# Compute tests
knueppel_test(u1)
berkowitz_test(u1)

# Simulate data that violates the null
u2 <- pnorm(arima.sim(model = list(ar = .5), n = 2000, sd = 1))
# Show PIT histogram
hist(u2)
# Compute tests
knueppel_test(u2)
berkowitz_test(u2)
}
\references{
Knueppel, M. (2015): `Evaluating the Calibration of Multi-Step-Ahead Density Forecasts Using Raw Moments'.
Journal of Business and Economic Statistics 33, 270-281. \doi{10.1080/07350015.2014.948175}
}
\seealso{
See \code{berkowitz_test for another test that applies to the same null hypothesis and data setup.}
}
