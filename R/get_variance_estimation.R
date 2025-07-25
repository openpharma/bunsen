#' Calculate the variance (SE) of the marginal treatment effects (hazard ratio) adjusting covariates in clinical trials
#'
#' @description
#' Estimate the standard error or variance of the marginal treatment effects using nonparametric bootstrap. Currently,
#' this only supports clustermq for parallel computation.
#' @details
#' If clustermq is not available, we suggest building your own bootstrap like boot and doParallel by using the function -- \link[bunsen]{get_point_estimate}.
#' This can also get you the SE or variance estimates. If you only run this function, you need to have cox_censor and cox_event in the environment.
#' @param cox_event Object. A coxph model using the survival time and survival status.
#' @param cox_censor Object. A coxph model using the survival time and 1-survival status.
#' @param trt Character. Variable name of the treatment assignment. Only support two arm trial at the moment.
#' @param data A data frame used for cox_event and cox_censor.
#' @param n.boot Numeric. Number of bootstrap.
#' @param M Numeric. The number of simulated counterfactual patients. Suggest to set above 1,000,000 to get robust estimation but it is time comsuming,
#' @param seed Numeric. Random seed for simulation.
#' @param cpp Bool. True for using C++ optimization. False for not using C++ optimization. This requires cpp package installed.
#' @param control Named list. A list containing control parameters, including memory of remote workers, whether to use nested parallel computation or local multiprocess, number of remote workers/jobs, etc. See details of \link[bunsen]{clmqControl}.
#' @param verbose Bool. Print status messages. Default: TRUE
#' @return A vector containing SE and 95% CI.
#'
#' @references
#' Daniel R, Zhang J, Farewell D. Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors
#' after adjustment for different covariate sets.
#' Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#' @export
#' @importFrom survival basehaz
#' @importFrom clustermq Q
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.formula coef model.frame na.omit predict quantile rbinom sd update var
#' @importFrom utils data
#' @examples
#' \dontrun{
#' #Don't run as it requires LSF scheduler
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' cox_censor <- coxph(Surv(OS, 1 - os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' get_variance_estimation(cox_event, cox_censor,
#'   trt = "trt", data = oak,
#'   M = 1000, n.boot = 10, control = clmqControl(), cpp = FALSE
#' )
#' }
get_variance_estimation <- function(cox_event, cox_censor, trt, data, M, n.boot, seed = NULL, cpp = TRUE,
                                    control = clmqControl(),verbose=TRUE) {
  if (is.null(seed)) seed <- Sys.time()
  if(verbose) cat(paste0("Calculating SE in clustermq using bootstrap N = ", n.boot, "...\n"))

  options(clustermq.scheduler = "LSF")


  out <- Q(
    fun = .fx.clsmq.cpp, i = 1:n.boot, n_jobs = control$n_jobs,
    const = list(
      memory = control$memory, data = data,
      cox_event = cox_event, cox_censor = cox_censor,
      trt = trt, local = control$local_se, clmq = control$clmq_se,
      M = M, cpp = cpp, seed = seed
    ),
    export = list(
      get_point_estimate = get_point_estimate,
      calculate_statistics = calculate_statistics,
      calculate_trt_effect = calculate_trt_effect,
      simulate_counterfactuals = simulate_counterfactuals,
      clmqControl = clmqControl,
      .fx_clsmq_simcoun = .fx_clsmq_simcoun
    ),
    template = list(cores = control$local_cores), pkgs = c("survival", "Rcpp", "clustermq"),verbose=verbose
  )
  hr_se <- do.call(c, out)
  se <- sd(hr_se)
  lb <- quantile(hr_se, prob = .025)
  ub <- quantile(hr_se, prob = .975)
  return(output = c(se = se, lb, ub))
}


.fx.clsmq.cpp <- function(i, data = data, cox_event = cox_event, cox_censor = cox_censor, trt = trt, M = M, seed = seed, cpp = cpp,
                          memory = memory, local = local, clmq = clmq) {
  tmp <- sample(1:nrow(data), size = nrow(data), replace = TRUE)

  tmp.dt <- data[tmp, ]

  cox_event_tmp <- update(cox_event, data = tmp.dt)

  cox_censor_tmp <- update(cox_censor, data = tmp.dt)

  hr_tmp <- get_point_estimate(
    trt = trt, cox_event = cox_event_tmp, cox_censor = cox_censor_tmp, data = tmp.dt, M = M, seed = seed, cpp = cpp,
    control = clmqControl(clmq_hr = clmq, clmq_local = local, memory = memory)

  )
  return(hr_tmp)
}
