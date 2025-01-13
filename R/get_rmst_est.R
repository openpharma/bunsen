#'Calculate the marginal restricted mean survival time (RMST) when adjusting covariates in clinical trials
#'
#'@description
#'Estimate the marginal RMST (point estimate) using the Karrison et al.(2018). Standard errors (SE) were estimated using the methods from Zucker (1998),
#'Chen and Tsiatis (2001), and Wei et al.(2023). We implemented both nonparametric(bootstrap) and parametric methods(delta) for SE.
#'
#'@details
#'Restricted mean survival time is a measure of average survival time up to a specified time point. We adopted the methods from Karrison et al.(2018) for
#'estimating the marginal RMST when adjusting covariates. For the SE, both nonparametric bootstrap and delta method are good for estimation.
#'We used the delta estimation from Zucker (1998) but we also included an additional variance component which treats the covariates as random as described in Chen and Tsiatis (2001).
#'
#'@param time A vector containing the event time of the sample.
#'@param status A vector containing the survival status of the sample.
#'@param arm A vector indicating the treatment assignment. 1 for treatment group. 0 for placebo group.
#'@param covariates A data frame containing the covariates. If covariates is NULL, unadjusted RMST is returned.
#'@param tau Numeric. A value for the restricted time or the pre-specified cutoff time point.
#'@param SE Character. If SE = 'boot', SE was estimated using nonparametric bootstrap. If 'delta', SE was estimated using the delta method. Default is 'delta'.
#'@param n.boot Numeric. Number of bootstrap used. Only used if SE = 'boot'.
#'
#'@return A list including marginal RMST and SE.
#'
#'@references
#'-  Karrison T, Kocherginsky M. Restricted mean survival time: Does covariate adjustment improve precision in randomized clinical trials? Clinical Trials. 2018;15(2):178-188. doi:10.1177/1740774518759281
#'-  Zucker, D. M. (1998). Restricted Mean Life with Covariates: Modification and Extension of a Useful Survival Analysis Method. Journal of the American Statistical Association, 93(442), 702–709. https://doi.org/10.1080/01621459.1998.10473722
#'-  Wei, J., Xu, J., Bornkamp, B., Lin, R., Tian, H., Xi, D., … Roychoudhury, S. (2024). Conditional and Unconditional Treatment Effects in Randomized Clinical Trials: Estimands, Estimation, and Interpretation. Statistics in Biopharmaceutical Research, 16(3), 371–381. https://doi.org/10.1080/19466315.2023.2292774
#'-  Chen, P. and Tsiatis, A. (2001), “Causal Inference on the Difference of the Restricted Mean Lifetime Between Two Groups,” Biometrics; 57: 1030–1038. DOI: 10.1111/j.0006-341x.2001.01030.x.
#'@export
#'
#'@importFrom survival coxph Surv strata
#'@importFrom boot boot
#'@examples
#'
#'data('oak')
#'tau=26
#'time=oak$OS
#'status=oak$os.status
#'arm=oak$trt
#'covariates=oak[,c("btmb","pdl1")]
#'get_rmst_est(time, status, arm, covariates,tau,SE='delta')
#'


get_rmst_est=function(time, status, arm, covariates=NULL,tau,SE='delta',n.boot=1000){

  if(is.null(covariates)){

    results=rmst_unadjust(time,status,arm,tau)

  }else{
    covariates=as.data.frame(covariates)

    dt=as.data.frame(cbind(time,status, arm, covariates))

    fit=coxph(as.formula(paste0('Surv(time, status) ~ ',paste(names(covariates),collapse = '+'),'+strata(arm)')),data = dt)

    delta=rmst_point_estimate(fit,dt=dt,tau)

    output=delta$output

    if(SE=='boot'){

      out=boot(data=dt,.rmst_boot_fx,R=n.boot,fit=fit,tau=tau,covariates=covariates)

      se=sd(out$t,na.rm = T)
      lb=quantile(out$t, prob=.025,na.rm=T)
      ub=quantile(out$t, prob=.975,na.rm=T)
      out=c(SE=se,lb,ub)

    }else if(SE=='delta') {

      out=rmst_delta(fit,time, arm, covariates,tau,surv0=delta$surv0,surv1=delta$surv1,cumhaz0=delta$cumhaz0,cumhaz1=delta$cumhaz1)
      out=c(SE=out,`2.5%`=output-1.96*out,`97.5%`=output+1.96*out)

    }

    results=list(RMST=output,SE=out)
  }

  return(results)

}

.rmst_boot_fx=function(data,idx,fit,tau,covariates){

  fit_tmp <- update(fit,data=data[idx,])
  tmax=basehaz(fit_tmp)
  tmax0=max(tmax$time[tmax$strata==unique(tmax$strata)[1]])
  tmax1=max(tmax$time[tmax$strata==unique(tmax$strata)[2]])
  if(tmax0<tau |tmax1<tau) {
    output=NA
  }else{
    out=rmst_point_estimate(fit=fit_tmp,dt=data[idx,],tau=tau)
    output=out$output
  }
  return(output)
}
