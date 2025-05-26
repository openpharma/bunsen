#' Control of marginal HR estimation via clustermq
#'
#'
#' @description
#' Construct control structures for marginal HR estimation. Specifically, this is the control for parallel computation via clustermq. It provides the nested parallel computation via LSF (remote parallel computation within remote parallel computation) and local multiprocess within remote parallel computation.
#'
#' @details
#' The control function provides options to set the memory of each remote node, number of cpus used by each remote node, and computation approach (whether or not use the nested parallel computation or local multiprocess with parallel computation).
#' @param memory Numeric. Memory allocation for the remote workers.
#' @param local_se Bool. True for calculating SE using local multiprocess in remote workers. This is only useful when clmq_se = TRUE.
#' @param clmq_se Bool. True for using parallel computation (nested or local) via clustermq. This can be combined with local_se to calculate SE with nested parallel computation or local multiprocess. Nested parallel computation means double parallel computations -- each worker will do a parallel computation for \link[bunsen]{simulate_counterfactuals}. False for calculating SE only in remote workers without nested parallel computation and local multiprocess.
#' @param clmq_hr Bool. True for calculating point estimate (marginal HR) using parallel computation via clustermq.
#' @param clmq_local Bool. True for calculating point estimate (marginal HR) using local multiprocess in remote workers. This is only useful when clmq_hr = TRUE.
#' @param n_jobs Numeric. Number of remote workers via clustermq.
#' @param local_cores Numeric. Number of cores or processes used in local multiprocess. This is only useful when local_se or clmq_local = TRUE.
#'
#' @return A list containing the control arguments.
#' @export
#'
#' @examples
#' \donttest{
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' cox_censor <- coxph(Surv(OS, 1 - os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' get_marginal_effect(
#'   trt = "trt", cox_event = cox_event, cox_censor = cox_censor, SE = TRUE,
#'   M = 1000, n.boot = 10, data = oak, seed = 1, cpp = FALSE, control = clmqControl()
#' )
#' #
#' }
#'
clmqControl <- function(memory = 1024 * 32, local_se = FALSE, clmq_se = FALSE, clmq_hr = TRUE, clmq_local = FALSE, n_jobs = 100, local_cores = 1) {
  return(structure(
    list(
      memory = memory, local_se = local_se, clmq_se = clmq_se,
      clmq_hr = clmq_hr, clmq_local = clmq_local, n_jobs = n_jobs, local_cores = local_cores
    ),
    class = c("clmqControl")
  ))
}
