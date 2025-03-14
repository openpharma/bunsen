#' Calculate the potential outcomes using marginal survival causual curves
#'
#' @description
#' Using the marginal survival causal curves from \link[bunsen]{calculate_statistics} to simulate the potential outcomes.
#'
#' @details
#' The potential outcomes were simulated by using a Bernoulli distribution from rbinom() and marginal survival causal curves. If M is quite large, we
#' suggest to use C++ optimization to speed up the calculation.
#'
#' @param bh A data frame from the \link[survival]{basehaz} function. This is the baseline hazard for the coxph model.
#' @param surv_cond A vector containing the marginal causal survival curves from \link[bunsen]{calculate_statistics}.
#' Each number is the probability of the surviving the time window (t1,t2],... conditional on surviving the prior corresponding window.
#' @param M Numeric. The number of simulated counterfactual patients. Suggest to set above 1,000,000 to get robust estimation but it is time comsuming,
#' @param cpp Bool. True for using C++ optimization. False for not using C++ optimization. This requires cpp package installed.
#' @param loadcpp Bool. True for loading C++ optimization functions. Default is TRUE. This is only used when cpp = TRUE.

#' @return A list containing the simulated event time and simulated status indicator.
#'
#' @references
#' Daniel R, Zhang J, Farewell D. Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors
#' after adjustment for different covariate sets.
#' Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#' @export
#' @importFrom Rcpp sourceCpp
#' @examples
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' cox_censor <- coxph(Surv(OS, 1 - os.status) ~ trt + btmb + pdl1, data = oak)
#'
#' bh <- basehaz(cox_event, centered = FALSE)
#' s_condi <- calculate_statistics(model = cox_event, trt = "trt")
#' sim_out_1d <- simulate_counterfactuals(bh = bh, surv_cond = s_condi$surv_cond0, cpp = FALSE, M = 1000)
simulate_counterfactuals <- function(bh, surv_cond, M, cpp, loadcpp = TRUE) {
  if (cpp) {
    if (loadcpp) sourceCpp("./src/cpp_functions.cpp")
    newd <- rbinom_matrix_vec(nrows = nrow(bh), ncols = M, surv_cond)
    index_event_d <- firstZeroIndex(newd)
  } else {
    newd <- matrix(stats::rbinom(n = nrow(bh) * M, size = 1, prob = surv_cond), nrow = nrow(bh), ncol = M)
    index_event_d <- apply(newd, 2, function(x) which(x == 0)[1])
  }

  eventtime_d <- bh$time[index_event_d]

  eventtime_d[which(is.na(index_event_d))] <- max(bh$time)

  status_d <- rep(1, M)

  status_d[which(is.na(index_event_d))] <- 0
  return(list(eventtime = eventtime_d, status = status_d))
}
