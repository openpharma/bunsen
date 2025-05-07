#' Calculate the point estimate of the marginal restricted mean survival time (RMST) when adjusting covariates in clinical trials
#'
#' @description
#' Estimate the marginal RMST (point estimate) using the Karrison et al.(2018).
#'
#' @details
#' Restricted mean survival time is a measure of average survival time up to a specified time point. We adopted the methods from Karrison et al.(2018) for
#' estimating the marginal RMST when adjusting covariates.
#' @param fit A \link[survival]{coxph} object with strata(trt) in the model. See example.
#' @param dt A data frame used for the fit - coxph model including survival time, OS status, trt, and covariates.
#' @param tau Numeric. A value for the restricted time or the pre-specified cutoff time point.

#' @return A list containing the RMST, cumulative survival function, and cumulative hazard function.
#' \describe{
#'   \item{output}{Marginal RMST}
#'   \item{surv0}{Cumulative survival function for the placebo group}
#'   \item{cumhaz0}{Cumulative hazard function for the placebo group}
#'   \item{surv1}{Cumulative survival function for the treatment group}
#'   \item{cumhaz1}{Cumulative hazard function for the treatment group}
#' }
#' @references
#' -  Karrison T, Kocherginsky M. Restricted mean survival time: Does covariate adjustment improve precision in randomized clinical trials? Clinical Trials. 2018;15(2):178-188. doi:10.1177/1740774518759281
#' -  Zucker, D. M. (1998). Restricted Mean Life with Covariates: Modification and Extension of a Useful Survival Analysis Method. Journal of the American Statistical Association, 93(442), 702–709. https://doi.org/10.1080/01621459.1998.10473722
#' -  Wei, J., Xu, J., Bornkamp, B., Lin, R., Tian, H., Xi, D., … Roychoudhury, S. (2024). Conditional and Unconditional Treatment Effects in Randomized Clinical Trials: Estimands, Estimation, and Interpretation. Statistics in Biopharmaceutical Research, 16(3), 371–381. https://doi.org/10.1080/19466315.2023.2292774
#' -  Chen, P. and Tsiatis, A. (2001), “Causal Inference on the Difference of the Restricted Mean Lifetime Between Two Groups,” Biometrics; 57: 1030–1038. DOI: 10.1111/j.0006-341x.2001.01030.x.
#' @export
#'
#' @importFrom survival basehaz
#' @examples
#' library(survival)
#' data("oak")
#'
#' tau <- 26
#' time <- oak$OS
#' status <- oak$os.status
#' trt <- oak$trt
#' covariates <- oak[, c("btmb", "pdl1")]
#' dt <- as.data.frame(cbind(time, status, trt, covariates))
#' fit <- coxph(Surv(time, status) ~ btmb + pdl1 + strata(trt), data = dt)
#' delta <- rmst_point_estimate(fit, dt = dt, tau)
#' delta$output

rmst_point_estimate <- function(fit, dt, tau) {
  cumhaz <- basehaz(fit, centered = FALSE)

  grp <- unique(cumhaz$strata)

  cumhaz0 <- cumhaz[cumhaz$strata == grp[1] & cumhaz$time <= tau, ]

  cumhaz1 <- cumhaz[cumhaz$strata == grp[2] & cumhaz$time <= tau, ]

  pred <- predict(fit, newdata = dt, type = "risk", reference = "zero")

  surv0 <- exp(-pred %*% t(cumhaz0$hazard))
  surv0_mean <- apply(surv0, 2, mean)

  d_t0 <- diff(c(0, cumhaz0$time, tau))

  mu0 <- sum(c(1, surv0_mean) * d_t0)

  surv1 <- exp(-pred %*% t(cumhaz1$hazard))
  surv1_mean <- apply(surv1, 2, mean)

  d_t1 <- diff(c(0, cumhaz1$time, tau))

  mu1 <- sum(c(1, surv1_mean) * d_t1)

  output <- mu1 - mu0

  return(list(
    output = output,
    surv0 = surv0,
    cumhaz0 = cumhaz0,
    surv1 = surv1,
    cumhaz1 = cumhaz1
  ))

}
