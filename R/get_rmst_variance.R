#'Calculate the variance of the marginal restricted mean survival time (RMST) when adjusting covariates using the delta method
#'
#'@description
#'Standard errors (SE) were estimated using the delta methods from Zucker (1998),
#'Chen and Tsiatis (2001), and Wei et al.(2023).
#'
#'@details
#'Restricted mean survival time is a measure of average survival time up to a specified time point. We adopted the methods from Karrison et al.(2018) for
#'estimating the marginal RMST when adjusting covariates. For the SE, both nonparametric bootstrap and delta method are good for estimation.
#'We used the delta estimation from Zucker (1998) but we also included an additional variance component which treats the covariates as random as described in Chen and Tsiatis (2001).
#'@param fit A \link[survival]{coxph} object with strata(arm) in the model. See example.
#'@param time A vector containing the event time of the sample.
#'@param arm A vector indicating the treatment assignment. 1 for treatment group. 0 for placebo group.
#'@param covariates A data frame containing the covariates.
#'@param tau Numeric. A value for the restricted time or the pre-specified cutoff time point.
#'@param surv0 A vector containing the cumulative survival function for the placebo group or arm0.
#'@param surv1 A vector containing the cumulative survival function for the treatment group or arm1.
#'@param cumhaz0 A data frame containing the cumulative hazard function for the placebo group or arm0.
#'@param cumhaz1 A data frame containing the cumulative hazard function for the placebo group or arm1.
#'
#'@return A value of the SE.
#'
#'@references
#'-  Karrison T, Kocherginsky M. Restricted mean survival time: Does covariate adjustment improve precision in randomized clinical trials? Clinical Trials. 2018;15(2):178-188. doi:10.1177/1740774518759281
#'-  Zucker, D. M. (1998). Restricted Mean Life with Covariates: Modification and Extension of a Useful Survival Analysis Method. Journal of the American Statistical Association, 93(442), 702–709. https://doi.org/10.1080/01621459.1998.10473722
#'-  Wei, J., Xu, J., Bornkamp, B., Lin, R., Tian, H., Xi, D., … Roychoudhury, S. (2024). Conditional and Unconditional Treatment Effects in Randomized Clinical Trials: Estimands, Estimation, and Interpretation. Statistics in Biopharmaceutical Research, 16(3), 371–381. https://doi.org/10.1080/19466315.2023.2292774
#'-  Chen, P. and Tsiatis, A. (2001), “Causal Inference on the Difference of the Restricted Mean Lifetime Between Two Groups,” Biometrics; 57: 1030–1038. DOI: 10.1111/j.0006-341x.2001.01030.x.
#'@export

#'@examples
#'library(survival)
#'data('oak')
#'
#'tau=26
#'time=oak$OS
#'status=oak$os.status
#'arm=oak$trt
#'covariates=oak[,c("btmb","pdl1")]
#'
#'dt=as.data.frame(cbind(time,status, arm, covariates))
#'fit=coxph(Surv(time, status) ~ btmb+pdl1+strata(arm),data = dt)
#'delta=rmst_point_estimate(fit,dt=dt,tau)
#'rmst_delta(fit,time, arm, covariates,tau,surv0=delta$surv0,surv1=delta$surv1,cumhaz0=delta$cumhaz0,cumhaz1=delta$cumhaz1)



rmst_delta=function(fit,time, arm, covariates,tau,surv0,surv1,cumhaz0,cumhaz1){
  sigma=fit$var*fit$n

  t0 = time[arm==0]

  t1 = time[arm==1]

  n0 = length(t0)

  n1 = length(t1)

  x0 = covariates[arm==0,]

  x1 = covariates[arm==1,]

  #arm 0
  s0_c = rep(0, length(cumhaz0$time))

  s1_c = matrix(0, nrow=length(cumhaz0$time), ncol=ncol(x0))

  zb = predict(fit,type='lp',reference = 'zero')

  zb0=zb[arm==0]

  zb1=zb[arm==1]

  for(i in 1:length(cumhaz0$time)){


    s0_c[i] <- sum((t0>=cumhaz0$time[i])*exp(zb0))/n0

    s1_c[i, ] <- colSums((t0>=cumhaz0$time[i])*exp(zb0)*x0)/n0

  }

  ht_c <- apply(s1_c/s0_c^2, 2, cumsum)/n0

  vt_c <- cumsum(1/s0_c^2)/n0


  gt_c <- apply(exp(zb)*surv0, 2, mean)


  pi_c <- matrix(0, nrow=length(cumhaz0$time), ncol=ncol(covariates))

  for(i in 1:ncol(covariates)){

    pi_c[, i] <- apply(covariates[, i]*exp(zb)*surv0, 2, mean)

  }



  mat_dt <- diff(c(cumhaz0$time, tau)) %*% t(diff(c(cumhaz0$time, tau)))

  mat_gt_c <- gt_c %*% t(gt_c)

  mat_vt_c <- matrix(0, nrow=length(cumhaz0$time), ncol=length(cumhaz0$time))

  for(i in 1:nrow(mat_vt_c)){

    for(j in 1:nrow(mat_vt_c)){

      id <- min(i, j)

      mat_vt_c[i, j] <- vt_c[id]

    }

  }

  omga_c <- sum(mat_vt_c * mat_gt_c * mat_dt)


  phi_c <- colSums((gt_c*ht_c - cumsum(diff(c(0,cumhaz0$hazard)))*pi_c)*diff(c(cumhaz0$time, tau)))


  var_c1 <- as.numeric(omga_c/n0 + phi_c %*% sigma %*% phi_c/fit$n)


  # variance induced by covariates

  AUC_c <- rep(NA, fit$n)

  for(i in 1:fit$n){

    AUC_c[i] <- sum(c(1, surv0[i, ])*diff(c(0, cumhaz0$time, tau)))

  }

  var_c2 <- mean((AUC_c-mean(AUC_c))^2)


  #arm 1
  s0_t = rep(0, length(cumhaz1$time))

  s1_t = matrix(0, nrow=length(cumhaz1$time), ncol=ncol(x1))


  for(i in 1:length(cumhaz1$time)){


    s0_t[i] <- sum((t1>=cumhaz1$time[i])*exp(zb1))/n1

    s1_t[i, ] <- colSums((t1>=cumhaz1$time[i])*exp(zb1)*x1)/n1

  }

  ht_t <- apply(s1_t/s0_t^2, 2, cumsum)/n1

  vt_t <- cumsum(1/s0_t^2)/n1


  gt_t <- apply(exp(zb)*surv1, 2, mean)


  pi_t <- matrix(0, nrow=length(cumhaz1$time), ncol=ncol(covariates))

  for(i in 1:ncol(covariates)){

    pi_t[, i] <- apply(covariates[, i]*exp(zb)*surv1, 2, mean)

  }



  mat_dt_t <- diff(c(cumhaz1$time, tau)) %*% t(diff(c(cumhaz1$time, tau)))

  mat_gt_t <- gt_t %*% t(gt_t)

  mat_vt_t <- matrix(0, nrow=length(cumhaz1$time), ncol=length(cumhaz1$time))

  for(i in 1:nrow(mat_vt_t)){

    for(j in 1:nrow(mat_vt_t)){

      id <- min(i, j)

      mat_vt_t[i, j] <- vt_t[id]

    }

  }

  omga_t <- sum(mat_vt_t * mat_gt_t * mat_dt_t)


  phi_t <- colSums((gt_t*ht_t - cumsum(diff(c(0,cumhaz1$hazard)))*pi_t)*diff(c(cumhaz1$time, tau)))


  var_t1 <- as.numeric(omga_t/n1 + phi_t %*% sigma %*% phi_t/fit$n)


  # variance induced by covariates

  AUC_t <- rep(NA, fit$n)

  for(i in 1:fit$n){

    AUC_t[i] <- sum(c(1, surv1[i, ])*diff(c(0, cumhaz1$time, tau)))

  }

  var_t2 <- mean((AUC_t-mean(AUC_t))^2)

  AUC_delta <- AUC_t-AUC_c

  cov_tc <- (fit$n-1)/fit$n*var(AUC_delta)


  se_c <- sqrt(var_c1 + var_c2/fit$n)

  se_t <- sqrt(var_t1 + var_t2/fit$n)


  var_d <- as.numeric(omga_c/n0 + omga_t/n1 + (phi_t-phi_c)%*%sigma%*%(phi_t-phi_c)/fit$n + cov_tc/fit$n)


  return(sqrt(var_d))


}


