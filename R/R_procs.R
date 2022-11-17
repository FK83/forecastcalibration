#' @useDynLib forecastcalibration, .registration = TRUE
#' @importFrom stats dnorm lm optim optimize pchisq qnorm dist pnorm
#' @importFrom Rcpp evalCpp
#' @importFrom scoringRules es_sample

# log likelihood for unconstrained AR(1) model
ll_ar1_unconstrained <- function(y, theta){
  # intercept parameter
  nu <- theta[1]
  # autoregressive parameter (between -1 and 1)
  a <- tanh(theta[2])
  # conditional standard deviation (aka residual standard deviation, > 0)
  s <- exp(theta[3])
  # vector of mean forecasts (t = 2, ..., T)
  m <- nu + a*y[-length(y)]
  # unconditional mean and variance (t = 1)
  um <- nu/(1-a)
  uv <- (s^2)/(1-a^2)
  # neg log likelihood for t = 1
  l1 <- -dnorm(y[1], mean = um, sd = sqrt(uv), log = TRUE)
  # neg log likelihood for = 2, ..., T
  l2 <- -sum(dnorm(y[-1], mean = m, sd = s, log = TRUE))
  # return complete neg log likelihood
  return(l1 + l2)
}
# log likelihood for constrained AR(1) model
ll_ar1_constrained <- function(y, a1){
  # conditional standard deviation
  s <- sqrt(1-a1^2)
  # conditional mean
  m <- a1*y[-length(y)]
  # unconditional mean and variance (for t = 1) given by 0 and 1
  l1 <- -dnorm(y[1], mean = 0, sd = 1, log = TRUE)
  # conditional mean and standard deviation (t = 2, ..., T)
  l2 <- -sum(dnorm(y[-1], mean = m, sd = s, log = TRUE))
  # return complete neg log lik
  return(l1 + l2)
}

#' Berkowitz-style test for calibration of multi step forecast distributions
#' @param u vector of PIT values
#' @return p-value for the null hypothesis of calibration
#' @details Under the null of calibration, the unconditional distribution of u is standard uniform.
#' However, the elements of u may be serially correlated under the null. The present variant of the Berkowitz
#' test is described in Knueppel (2015). As described there, the test is not exact since the AR(1) model under the alternative
#' may be misspecified. However, the test often works well in practice.
#' @references
#' Berkowitz, J. (2001): `Testing Density Forecasts, With Applications to Risk Management'.
#' Journal of Business and Economic Statistics 19, 465-474. \doi{10.1198/07350010152596718}
#'
#' Knueppel, M. (2015): `Evaluating the Calibration of Multi-Step-Ahead Density Forecasts Using Raw Moments'.
#' Journal of Business and Economic Statistics 33, 270-281. \doi{10.1080/07350015.2014.948175}
#'
#' @export
berkowitz_test <- function(u){

  if (any((u == 0) | (u == 1))){

    stop("berkowitz_test cannot handle PIT-values of exactly zero or one. knueppel_test is more appropriate in this case.")

  }

  y <- qnorm(u)

  # fit unconstrained AR(1) model
  opt_u <- optim(c(0, atanh(.7), 0), fn = ll_ar1_unconstrained, y = y)
  # likelihood value
  ll_u <- -opt_u$value

  # fit constrained AR(1) model
  opt_c <- optimize(f = ll_ar1_constrained, y = y, interval = c(-1, 1))
  ll_c <- -opt_c$objective

  # likelihood ratio test
  lr_stat <- 2*(ll_u - ll_c)
  lr_p <- 1-pchisq(lr_stat, df = 2)

  list(lr_stat = lr_stat, lr_p = lr_p)

}

# Function to compute HAC covariance matrix for (possibly) multivariate
# response u. Translation of Matlab code by Malte Knueppel.
vHAC <- function(u, prewhite = FALSE, k = NULL, meth = "qs"){

  if (!is.matrix(u)) u <- as.matrix(u)

  n <- nrow(u)
  nreg <- ncol(u)
  rho <- sigma <- matrix(0, nreg, 1)

  if (prewhite == TRUE){
    # do VAR(1) prewhitening
    reg.x <- matrix(u[1:(n-1), ], n-1, nreg)
    reg.y <- matrix(u[2:n, ], n-1, nreg)
    aux <- lm(reg.y~reg.x-1, )
    beta <- matrix(unname(aux$coefficients), nreg, nreg)
    v <- matrix(unname(aux$residuals), n-1, nreg)
  } else {
    v <- u
    beta <- matrix(0, nreg, nreg)
  }
  nv <- nrow(v)

  # choose nr of lags (if not provided)
  if (is.null(k)){

    for (i in 1:nreg){
      aux <- lm(v[2:nv, i]~v[1:(nv-1), i]-1)
      rho[i] <- unname(aux$coefficients)
      sigma[i] <- sum(unname(aux$residuals)^2) / nv
    }

    if (meth == "qs"){
      top <- sum( (4*(rho^2) * (sigma^2)) / ((1-rho)^8) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- ceiling(1.3221*((top/bot)*n)^(0.2))
    } else {
      top <- sum( (4*(rho^2) * (sigma^2)) / (((1-rho)^6)*((1+rho)^2)) )
      bot <- sum( (sigma^2) / ((1-rho)^4) )
      k <- min(c(ceiling(1.1447*((top/bot)*n)^(1/3)), round(0.5*n)))
    }

  }

  # compute HAC
  vcv <- (t(v) %*% v) / (n-1)

  if (k > 0){
    if (meth == "qs"){
      del <- ((6 * pi)/(5 * k)) * (1:(n-1))
      w <- 3 * (sin(del) / del - cos(del)) / (del^2)
      if (prewhite == FALSE){
        mlag <- n - 1
      } else {
        mlag <- nv - 1
      }
    } else {
      w <- 1 - (1:k)/(k+1)
      mlag <- k
    }
    for (i in 1:mlag){
      cov <- t(v[(i+1):nv, , drop = FALSE]) %*% (v[1:(nv-i), , drop = FALSE]) / (n-1)
      vcv <- vcv + w[i]*(cov + t(cov))
    }
  }

  d <- solve(diag(nreg) - t(beta))
  hac <- d %*% vcv %*% t(d)

  return(list(hac = hac, k = k))

}

#' Compute autocorrelation-robust t statistic
#'
#' @importFrom sandwich NeweyWest bwNeweyWest
#' @importFrom magrittr "%>%"
#' @export
#' @param x vector of data (with mean = 0 under null hypothesis)
#' @param variance_under_h0 logical, should null hypothesis be imposed for variance estimation?
#' @param implementation string, either "MK" or "sandwich". The former is the implementation of Knueppel (2015), the
#' latter uses the function \code{NeweyWest} from the sandwich package
#' @param ... additional parameters that can be passed to the robust variance function
#' @examples
#' # Simulate AR(1) model with zero mean
#' y <- arima.sim(model = list(ar = .5), n = 500)
#' # Different test variants
#' t_hac(y, implementation = 'MK')
#' t_hac(y, implementation = 'MK', variance_under_h0 = TRUE)
#' t_hac(y, implementation = 'sandwich')
t_hac <- function(x, implementation = "MK", variance_under_h0 = FALSE,
                  ...){
  if (implementation == "sandwich" & variance_under_h0){
    stop("Variance under H0 option not implemented for sandwich variant")
  }
  # Compute mean
  m <- mean(x)
  if (implementation == "MK"){
    if (variance_under_h0){
      v <- vHAC(x, ...)
    } else {
      v <- vHAC(x-m, ...)
    }
    k <- v$k
    t <- m/sqrt(v$hac/length(x))
  } else if (implementation == "sandwich"){
    fit <- lm(x~1)
    # Standard deviation
    s <- NeweyWest(fit, ...) %>% unname %>% sqrt
    # Bandwidth (aka lag length)
    k <- bwNeweyWest(fit, ...)
    # t stat
    t <- m/s
  }
  # Return t-stat, lag length and (two-sided) p-value
  return(list(t = t, k = k, pval = 2*pnorm(-abs(t))))
}

#' Knueppel (2015) test for calibration of multi-step forecast distribution
#'
#' @param x series of PIT values. Under the null of calibration, unconditional
#' distribution of \code{x} is standard uniform.
#' @param lags,prewhite,meth Parameters for computing the autocorrelation-robust variance
#' that enters the test statistic.
#' @param moments Choice of moments of uniform distribution to be considered. String, either \code{"1234"} (first four moments, the default) or \code{"12"} (first two moments).
#' @details The variants of the test are called "$alpha_{1234}^{0}$ and "$alpha_{12}^{0}$ in Knueppel (2015).
#' @references
#' Knueppel, M. (2015): `Evaluating the Calibration of Multi-Step-Ahead Density Forecasts Using Raw Moments'.
#' Journal of Business and Economic Statistics 33, 270-281. \doi{10.1080/07350015.2014.948175}
#' @examples
#' set.seed(1)
#'
#' # Simulate data that satisfies the null
#' u1 <- pnorm(arima.sim(model = list(ar = .5), n = 2000, sd = sqrt(1-.5^2)))
#' # Show PIT histogram
#' hist(u1)
#' # Compute tests
#' knueppel_test(u1)
#' berkowitz_test(u1)
#'
#' # Simulate data that violates the null
#' u2 <- pnorm(arima.sim(model = list(ar = .5), n = 2000, sd = 1))
#' # Show PIT histogram
#' hist(u2)
#' # Compute tests
#' knueppel_test(u2)
#' berkowitz_test(u2)
#' @seealso See \code{berkowitz_test for another test that applies to the same null hypothesis and data setup.}
#' @export
knueppel_test <- function(x, lags = NULL, prewhite = FALSE, meth = "qs",
                          moments = "1234"){

  s_pit <- sqrt(12) * (x - 0.5)
  if (moments == "12"){
    z_orig <- cbind(s_pit, (s_pit^2-1))
    z1 <- z_orig[, 1, drop = FALSE]
    z2 <- z_orig[, 2, drop = FALSE]
    d_f <- 2
  } else if (moments == "1234"){
    z_orig <- cbind(s_pit, (s_pit^2-1),  (s_pit^3), (s_pit^4-1.8))
    z1 <- z_orig[, c(1, 3)]
    z2 <- z_orig[, c(2, 4)]
    d_f <- 4
  }

  y1 <- matrix(colSums(z1)) * (nrow(z1)^(-0.5))
  tmp1 <- vHAC(z1, prewhite = prewhite, k = lags, meth = meth)
  phi1 <- tmp1$hac
  stat_odd <- t(y1) %*% solve(phi1) %*% y1

  y2 <- matrix(colSums(z2)) * (nrow(z2)^(-0.5))
  tmp2 <- vHAC(z2, prewhite = prewhite, k = lags, meth = meth)
  phi2 <- tmp2$hac
  stat_even <- t(y2) %*% solve(phi2) %*% y2

  stat <- stat_odd + stat_even
  pval <- 1 - pchisq(stat, d_f)

  return(list(stat = stat, pval = pval, nlags_odd = tmp1$k, nlags_even = tmp2$k,
              moments = moments))

}


#' Computations for score calibration test

#' @export
#' @param y realization, vector of length d (= dimension of forecast distribution)
#' @param x forecast sample, matrix with d rows (variables) and m columns (forecast draws).
#' @param exact should exact estimator be used? If \code{FALSE} (the default), split sample in two halfs.
#' @param randomize should sample split be randomized? Defaults to \code{FALSE}. Only relevant if \code{exact = FALSE}.
#' @param use_C should Cpp implementation be used? Defaults to \code{TRUE}; alternatively, use base R code based on \code{stats::dist}. Should yield the same results. C code is usually faster.
#' @param checks should inputs be checked for consistency? Defaults to \code{TRUE}.
#' @return list with three components, \code{u}, \code{d} and \code{es}. All are scalar numbers. \code{u} is the PIT value of the realized score within the simulated scores. \code{d} is the difference between the realized score and the expected score. \code{es} is the energy score.
#' @examples
#' set.seed(1)
#' # Draw five-dimensional forecast distribution and realization
#' x <- matrix(rnorm(5000), 5, 1000)
#' y <- rnorm(5)
#' # Compute score calibration quantities
#' score_calibration_es(y, x)
#' # Compare to scoringRules function (should yield the same
#' # energy score value if exact version is used)
#' scoringRules::es_sample(y, x)
score_calibration_es <- function(y, x, exact = FALSE, randomize = FALSE,
                                 use_C = TRUE, checks = TRUE){
  if (checks){
    input_checks_multiv(y, x)
  }
  if (exact){
    # Compute exact estimator
    if (use_C){
      aux <- esC_sample(y, x)
      aux_x <- aux$avgdist %>% as.numeric
      aux_y <- aux$eyx
    } else {
      aux_x <- x %>% t %>% dist %>% as.matrix %>% colMeans
      aux_y <- ((y-x)^2) %>% colSums %>% sqrt %>% mean
    }
  } else {
    if (!use_C){
      warning("Base R implementation available for exact estimator only. Input ignored.")
    }
    m <- ncol(x)
    m1 <- floor(.5*m)
    if (randomize){
      # randomize split samples
      sel_ind1 <- sample.int(m, size = m1, replace = FALSE)
      sel_ind2 <- setdiff(1:m, sel_ind1)
    } else {
      # split sample in two halfs
      sel_ind1 <- 1:m1
      sel_ind2 <- (m1+1):m
    }
    x1 <- x[,sel_ind1,drop = FALSE]
    x2 <- x[,sel_ind2,drop = FALSE]
    # run C code for computations
    aux <- esC_split(y, x1, x2)
    aux_x <- aux$avgdist %>% as.numeric
    aux_y <- aux$eyx
  }
  # PIT of actual score within n scores
  u <- mean(aux_x < aux_y)
  # diff. between actual and expected score
  d <- aux_y - mean(aux_x)
  # energy score
  es <- aux_y - .5*mean(aux_x)
  # return list
  list(u = u, d = d, es = es)
}

input_checks_multiv <- function(y, x){
  stopifnot(is.numeric(y))
  stopifnot(is.matrix(x))
  y_okay <- is.null(dim(y))
  if (!y_okay){
    stop("y must be a vector")
  }
  d <- length(y)
  if ( (nrow(x) != d) ){
    stop("forecast samples have wrong number of rows")
  }
}
