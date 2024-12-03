get_rmst_est=function(time, status, arm, covariates,tau,SE='delta',n.boot=1000){

  if(is.null(covariates)){

    results=.rmst_unadjust(time,status,arm,tau)

  }else{
    covariates=as.data.frame(covariates)

    data=as.data.frame(cbind(time,status, arm, covariates))

    fit=coxph(as.formula(paste0('Surv(time, status) ~ ',paste(names(covariates),collapse = '+'),'+strata(arm)')),data = data)

    delta=rmst_point_estimate(fit,tau)

    output=delta$output

    if(SE=='boot'){

      out=boot(data=data,rmst_boot_fx,R=n.boot,fit=fit,tau=tau)

      se=sd(out$t,na.rm = T)
      lb=quantile(out$t, prob=.025,na.rm=T)
      ub=quantile(out$t, prob=.975,na.rm=T)

      out=list(c(SE=se,lb,ub),t=out$t)

    }else if(SE=='delta') {

      out=rmst_delta(fit,time, arm, covariates,tau,surv0=delta$surv0,surv1=delta$surv1,cumhaz0=delta$cumhaz0,cumhaz1=delta$cumhaz1)

    }

    results=list(RMST=output,SE=out)
  }

  return(results)

}

rmst_point_estimate=function(fit,tau){

  cumhaz=basehaz(fit,centered = FALSE)

  grp <- unique(cumhaz$strata)

  cumhaz0 <- cumhaz[cumhaz$strata==grp[1] & cumhaz$time<=tau, ]

  cumhaz1 <- cumhaz[cumhaz$strata==grp[2] & cumhaz$time<=tau, ]

  pred=predict(fit,type='risk',reference = 'zero')

  surv0=exp(-pred%*%t(cumhaz0$hazard))
  surv0_mean= apply(surv0, 2, mean)

  d_t0=diff(c(0,cumhaz0$time,tau))

  mu0=sum(c(1,surv0_mean)*d_t0)

  surv1=exp(-pred%*%t(cumhaz1$hazard))
  surv1_mean=apply(surv1, 2, mean)

  d_t1=diff(c(0,cumhaz1$time,tau))

  mu1=sum(c(1,surv1_mean)*d_t1)

  output=mu1-mu0
  return(list(output=output,surv0=surv0,cumhaz0=cumhaz0,surv1=surv1,cumhaz1=cumhaz1))
}
