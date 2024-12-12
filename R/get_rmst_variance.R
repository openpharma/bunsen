

rmst_boot_fx=function(data,idx,fit,tau,covariates){

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


