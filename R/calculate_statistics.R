#' Estimate the marginal causal survival curves using potential outcome framework
#'
#' @description
#' Estimate the marginal causal survival curves for simulating time-to-event data in a discrete manner
#' based on the methods from Daniel et al.(2020).
#'
#'
#'
#' @details
#' If the study period for the original data is divided into discrete windows,
#' defined by the event times in the original data, at time t0 = 0, everyone in the simulated data is still a survivor.
#' S(x) is the estimated survival function. By the end of the window (0,t1], a proportion S(t1) still survives.
#' The conditional probability of surviving the next window, (t1,t2], conditional on surviving the first window,
#' is S(t2)∕S(t1), and so on. This function returns the S(t2)/S(t1) in series.
#'
#' @param model A fitted \link[survival]{coxph} model. This should be a coxph event model or censoring model.
#' @param trt Character. Name of the treatment assignment variable.

#' @return Two vectors containing the marginal causal survival curves for treatment arms (1 for treatment arm; 0 for control arm/placebo). Each number is the probability of the surviving the time window
#' (t1,t2],... conditional on surviving the prior corresponding window.
#'
#' @references
#' Daniel R, Zhang J, Farewell D. Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors
#' after adjustment for different covariate sets.
#' Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#' @export
#'
#' @examples
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' calculate_statistics(model = cox_event, trt = "trt")
#'
calculate_statistics <- function(model, trt) {
  bh <- basehaz(model, centered = FALSE)

  dt0 <- dt1 <- model.frame(model)

  ntrt <- sort(unlist(unique(dt0[trt])))

  dt0[trt] <- ntrt[1]

  dt1[trt] <- ntrt[2]

  pred <- predict(model, newdata = dt0, type = "risk", reference = "zero")
  survf <- exp(-pred %*% t(bh$hazard))

  survf_mean <- colMeans(survf)
  survf_1 <- c(1, survf_mean[1:(length(survf_mean) - 1)])
  surv_cond0 <- survf_mean / survf_1

  pred <- predict(model, newdata = dt1, type = "risk", reference = "zero")
  survf <- exp(-pred %*% t(bh$hazard))

  survf_mean <- colMeans(survf)
  survf_1 <- c(1, survf_mean[1:(length(survf_mean) - 1)])
  surv_cond1 <- survf_mean / survf_1

  output <- list(surv_cond0 = surv_cond0, surv_cond1 = surv_cond1)


  return(output)
}
