#' Calculate the marginal treatment effects for hazard ratio (HR) adjusting covariates in clinical trials
#'
#' @description
#' Use the simulation approach from Daniel et al. (2020) to estimate the marginal HR when adjusting covariates.
#' Standard error of marginal HR was estimated via the nonparametric bootstrap.
#' The standard error of HR will converge with the increase of bootstrap. We suggest setting a large M
#' for a more robust estimation. This function uses the parallel computation via remote LSF and local multiprocess. Additional
#' feature also includes the C++ optimization that can speed up the calculation.
#'
#' @details
#' In clinical trials, adjusting covariates like prognostic factors in the main analysis can increase the precision
#' resulting in smaller SE. However, adjusting covariates in nonlinear models changes the target estimands, e.g.
#' from marginal treatment effects to conditional treatment effects. This function has implemented the methods of Daniel et al. (2020)
#' to estimate the marginal treatment effects when adjusting covariates. Increasing the M - the number of bootstrap can significantly increase
#' the computation time. Hence, we introduced C++ optimization and parallel computation to speed up the calculation. It provides nested parallel
#' computation via LSF and remote parallel computation with local multiprocess.
#'
#' @param trt Character. Variable name of the treatment assignment. Only support two arm trial at the moment.
#' @param cox_event Object. A coxph model using the survival time and survival status.
#' @param cox_censor Object. A coxph model using the survival time and 1-survival status.
#' @param data A data frame used for cox_event and cox_censor.
#' @param M Numeric. The number of simulated counterfactual patients. Suggest to set a large number to get robust results, but this will be very time comsuming.
#' @param SE Bool. True for estimating SE. False for not estimating SE.
#' @param seed Numeric. Random seed for simulation.
#' @param cpp Bool. True for using C++ optimization. False for not using C++ optimization. This requires cpp package installed.
#' @param control Named list. A list containing control parameters, including memory of remote workers, whether to use nested parallel computation or local multiprocess, number of remote workers/jobs, etc. See details of \link[bunsen]{clmqControl}.
#' @param n.boot Numeric. Number of bootstrap used.

#'
#' @return A vector containing the marginal beta (logHR), standard error, and 95% CI.
#'
#' @references
#' Daniel R, Zhang J, Farewell D. Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors
#' after adjustment for different covariate sets.
#' Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#' @export
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' cox_censor <- coxph(Surv(OS, 1 - os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' get_marginal_effect(
#'   trt = "trt", cox_event = cox_event, cox_censor = cox_censor, SE = TRUE,
#'   M = 1000, n.boot = 10, data = oak, seed = 1, cpp = FALSE, control = clmqControl(n_jobs = 100)
#' )
#' #
#' }
get_marginal_effect <- function(trt, cox_event, cox_censor, data, M, SE = TRUE, seed = NULL, cpp = TRUE, n.boot = 1000,
                                control = clmqControl()) {
  sanitize_coxmodel(cox_event, trt)
  sanitize_coxmodel(cox_censor, trt)


  output <- get_point_estimate(
    trt = trt, cox_event = cox_event, cox_censor = cox_censor, data = data, M = M, seed = seed, cpp = cpp,
    control = control
  )

  if (SE) {
    SE <- get_variance_estimation(
      trt = trt, data = data, M = M, n.boot = n.boot, cpp = cpp, control = control,
      cox_event = cox_event, cox_censor = cox_censor, seed = seed
    )
    output <- c(beta = output, SE)
  }

  return(output)
}
