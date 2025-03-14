#' Calculate the unadjusted restricted mean survival time (RMST)
#'
#' @description
#' Estimate the unadjusted RMST (point estimate).

#' @param time A vector containing the event time of the sample.
#' @param status A vector containing the survival status of the sample.
#' @param arm A vector indicating the treatment assignment. 1 for treatment group. 0 for placebo group.
#' @param tau Numeric. A value for the restricted time or the pre-specified cutoff time point.

#' @return A data frame including the survival time for each arm and the difference. SE were also calculated.
#' \describe{
#'   \item{mu0}{Mean survival time for arm0}
#'   \item{se0}{SE of mu0}
#'   \item{mu1}{Mean survival time for arm1}
#'   \item{se1}{SE of mu1}
#'   \item{delta}{Difference between mu0 and mu1}
#'   \item{se_d}{SE of delta}
#' }
#' @export
#'
#' @examples
#'
#' data("oak")
#' tau <- 26
#' time <- oak$OS
#' status <- oak$os.status
#' arm <- oak$trt
#' covariates <- oak[, c("btmb", "pdl1")]
#' results <- rmst_unadjust(time, status, arm, tau)
#'
rmst_unadjust <- function(time, status, arm, tau) {
  # arm0

  t0 <- time[arm == 0]

  cn0 <- status[arm == 0]

  ft0 <- unique(sort(t0[cn0 == 1 & t0 <= tau]))

  nd0 <- sapply(ft0, function(x) sum(t0 == x & cn0 == 1))

  nrisk0 <- sapply(ft0, function(x) sum(t0 >= x))

  surv0 <- cumprod(1 - nd0 / nrisk0)


  # arm1

  t1 <- time[arm == 1]

  cn1 <- status[arm == 1]

  ft1 <- unique(sort(t1[cn1 == 1 & t1 <= tau]))

  nd1 <- sapply(ft1, function(x) sum(t1 == x & cn1 == 1))

  nrisk1 <- sapply(ft1, function(x) sum(t1 >= x))

  surv1 <- cumprod(1 - nd1 / nrisk1)


  # point estimate of RMST

  mu0 <- sum(c(1, surv0) * diff(c(0, ft0, tau)))

  mu1 <- sum(c(1, surv1) * diff(c(0, ft1, tau)))

  delta <- mu1 - mu0


  # estimate the variance

  V0i <- rep(NA, length(ft0))

  for (i in 1:length(ft0)) {
    if (nrisk0[i] == nd0[i]) {
      V0i[i] <- sum(surv0[ft0 >= ft0[i]] * diff(c(ft0[ft0 >= ft0[i]], tau)))^2 * nd0[i] / nrisk0[i]
    } else {
      V0i[i] <- sum(surv0[ft0 >= ft0[i]] * diff(c(ft0[ft0 >= ft0[i]], tau)))^2 * nd0[i] / nrisk0[i] / (nrisk0[i] - nd0[i])
    }
  }


  V1i <- rep(NA, length(ft1))

  for (i in 1:length(ft1)) {
    if (nrisk1[i] == nd1[i]) {
      V1i[i] <- sum(surv1[ft1 >= ft1[i]] * diff(c(ft1[ft1 >= ft1[i]], tau)))^2 * nd1[i] / nrisk1[i]
    } else {
      V1i[i] <- sum(surv1[ft1 >= ft1[i]] * diff(c(ft1[ft1 >= ft1[i]], tau)))^2 * nd1[i] / nrisk1[i] / (nrisk1[i] - nd1[i])
    }
  }


  V0 <- sum(V0i)

  V1 <- sum(V1i)



  data.frame(mu0, se0 = sqrt(V0), mu1, se1 = sqrt(V1), delta = delta, se_d = sqrt(V0 + V1))
}
