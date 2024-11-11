
get_rmst_estimate=function(fit,data,tau){

  cumhaz=basehaz(fit,centered = FALSE)

  grp <- unique(cumhaz$strata)

  cumhaz0 <- cumhaz[cumhaz$strata==grp[1] & cumhaz$time<=tau, ]

  cumhaz1 <- cumhaz[cumhaz$strata==grp[2] & cumhaz$time<=tau, ]

  pred=predict(fit,newdata = data,type='risk',reference = 'zero')

  surv0=exp(-pred%*%t(cumhaz0$hazard))
  surv0_mean= apply(surv0, 2, mean)

  d_t0=diff(c(0,cumhaz0$time,tau))

  mu0=sum(c(1,surv0_mean)*d_t0)

  surv1=exp(-pred%*%t(cumhaz1$hazard))
  surv1_mean=apply(surv1, 2, mean)

  d_t1=diff(c(0,cumhaz1$time,tau))

  mu1=sum(c(1,surv1_mean)*d_t1)

  delta=mu1-mu0

  return(delta)

}







