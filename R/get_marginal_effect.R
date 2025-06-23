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
#' @param verbose Bool. Print status messages. Default: TRUE
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
#' #Don't run as it requires LSF scheduler
#'
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
                                control = clmqControl(),verbose=TRUE) {
  sanitize_coxmodel(cox_event, trt)
  sanitize_coxmodel(cox_censor, trt)


  output <- get_point_estimate(
    trt = trt, cox_event = cox_event, cox_censor = cox_censor, data = data, M = M, seed = seed, cpp = cpp,
    control = control,verbose=verbose
  )

  ret <- list(
    trt=trt,
    cox_event=cox_event,
    cox_censor=cox_censor,
    data=data,
    M=M,
    seed=seed,
    cpp=cpp,
    control=control,
    beta=output
  )

  if (SE) {
    SE <- get_variance_estimation(
      trt = trt, data = data, M = M, n.boot = n.boot, cpp = cpp, control = control,
      cox_event = cox_event, cox_censor = cox_censor, seed = seed,verbose=verbose
    )
    ret$SE=SE
  }else{ret$SE=NA}

  class(ret) <- c('marginal_cox','bunsen')

  return(ret)
}


#' Summarizing the marginal treatment effects for hazard ratio (HR)
#'
#' summary method for class 'marginal_cox'.
#'
#' @param object an object of class 'marginal_cox'
#' @param ... Parameters for other methods.
#' @return No return value. This is called for its side effects.
#' @importFrom stats pnorm
#' @keywords internal
#' @export
summary.marginal_cox <- function(object,...){
  tmp=summary(object$cox_event)$coefficients
  e_cox=tmp[rownames(tmp)=='trt',]
  z_value <- ifelse(length(object$SE>1), as.numeric(object$beta / object$SE[1]),NA)
  p_value <- ifelse(length(object$SE>1),as.numeric(2 * (1 - pnorm(abs(z_value)))),NA)
  ret=structure(list(
    trt=object$trt,
    formula=object$cox_event$formula,
    original_cox=e_cox,
    M=object$M,
    N=nrow(object$data),
    nevent=object$cox_event$nevent,
    nevent_c=object$cox_censor$nevent,
    seed=object$seed,
    coef=object$beta,
    SE=object$SE,
    z_value=z_value,
    p_value=p_value

  ),class='summary.marginal_cox')
  return(ret)
}


#' Print the summary of marginal treatment effects for hazard ratio (HR)
#'
#' print method for class 'summary_marginal_cox'.
#'
#' @param x an object of class 'summary_marginal_cox'
#' @param ... Parameters for other methods.
#' @return No return value. This is called for its side effects.
#' @keywords internal
#' @export
print.summary.marginal_cox <- function(x,...){

  if(inherits(x,'summary.marginal_cox')){
    cat('Call:\n')
    print(x$formula)
    cat('Marginal treatment effect calculated by N =',x$M,'simulations\n')
    cat('Treatment variable:',x$trt,'------ Number of sample:',x$N,'\n')
    cat('Number of events in cox_event:',x$nevent,'\n')
    cat('Number of events in cox_censor:',x$nevent_c,'\n')
    cat('Random seed = ',x$seed,'\n')
    cat('Original treatment effect:\n')
    cat(sprintf("%-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",'','coef','exp(coef)','se(coef)','z','Pr(>|z|)'))
    cat(sprintf("%-10s  %-10f  %-10f  %-10f  %-10f  %-10f\n",
                x$trt,x$original_cox['coef'],x$original_cox['exp(coef)'],x$original_cox['se(coef)'],x$original_cox['z'],x$original_cox['Pr(>|z|)']))
    cat('Marginal treatment effect:\n')
    cat(sprintf("%-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",'','coef','exp(coef)','se(coef)','z','Pr(>|z|)'))
    cat(sprintf("%-10s  %-10f  %-10f  %-10f  %-10f  %-10f\n",x$trt,x$coef,exp(x$coef),x$SE[1],x$z_value,x$p_value))
    cat('95% CI of Marginal treatment effect (bootstrap):',sprintf('%.3f %s %.3f',x$SE[2],',',x$SE[3]))

  }

}


#' Print the marginal treatment effects for hazard ratio (HR)
#'
#' print method for class 'marginal_cox'.
#'
#' @param x an object of class 'marginal_cox'
#' @param ... Parameters for other methods.
#' @return No return value. This is called for its side effects.
#' @keywords internal
#' @export
print.marginal_cox <- function(x,...){

  if(inherits(x,'marginal_cox')){
    cat('Call:\n')
    print(x$cox_event$formula)
    cat('Marginal treatment effect calculated by N =',x$M,'simulations\n')
    cat('Number of sample:',nrow(x$data),'\n')
    cat(sprintf("%-10s  %-10s  %-10s  %-10s  %-10s  %-10s\n",'','coef','exp(coef)','se(coef)','2.5%','97.5%'))
    cat(sprintf("%-10s  %-10f  %-10f  %-10f  %-10f  %-10f\n",x$trt,x$beta,exp(x$beta),x$SE[1],x$SE[2],x$SE[3]))
    if(x$control$clmq_hr) cat('clustermq setting:\n','number of remote workers = ',x$control$n_jobs,', each worker has',x$control$local_cores,'core(s)\n')
    cat('Point estimate:',ifelse(x$control$clmq_hr,
                                 ifelse(x$control$clmq_local,
                                        'parallel computation (clustermq) with local multiprocess\n',
                                        'parallel computation (clustermq)\n'),
                                 'local environment\n'))
    if(length(x$SE)>1){
      cat('SE (bootstrap):',ifelse(x$control$clmq_se,
                                   ifelse(x$control$local_se,
                                          'parallel computation (clustermq) with local multiprocess',
                                          'nested parallel computation (clustermq)'),
                                   'parallel computation\n'))
      cat('95%CI estimated by bootstrap\n')
    }else{
      cat('SE: Not estimated\n')
    }


  }

}



