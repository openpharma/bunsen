.rmst_unadjust=function(time,status,arm,tau){

  # arm0

  t0 <- time[arm==0]

  cn0 <- status[arm==0]

  ft0 <- unique(sort(t0[cn0==1 & t0<=tau]))

  nd0 <- sapply(ft0, function(x) sum(t0==x & cn0==1))

  nrisk0 <- sapply(ft0, function(x) sum(t0>=x))

  surv0 <- cumprod(1-nd0/nrisk0)


  # arm1

  t1 <- time[arm==1]

  cn1 <- status[arm==1]

  ft1 <- unique(sort(t1[cn1==1 & t1<=tau]))

  nd1 <- sapply(ft1, function(x) sum(t1==x & cn1==1))

  nrisk1 <- sapply(ft1, function(x) sum(t1>=x))

  surv1 <- cumprod(1-nd1/nrisk1)


  # point estimate of RMST

  mu0 <- sum(c(1, surv0)*diff(c(0, ft0, tau)))

  mu1 <- sum(c(1, surv1)*diff(c(0, ft1, tau)))

  delta <- mu1-mu0


  # estimate the variance

  V0i <- rep(NA, length(ft0))

  for(i in 1:length(ft0)){

    if(nrisk0[i]==nd0[i]){
      V0i[i] <- sum(surv0[ft0>=ft0[i]]*diff(c(ft0[ft0>=ft0[i]], tau)))^2*nd0[i]/nrisk0[i]
    } else {
      V0i[i] <- sum(surv0[ft0>=ft0[i]]*diff(c(ft0[ft0>=ft0[i]], tau)))^2*nd0[i]/nrisk0[i]/(nrisk0[i]-nd0[i])
      }


  }


  V1i <- rep(NA, length(ft1))

  for(i in 1:length(ft1)){

    if(nrisk1[i]==nd1[i]) {
      V1i[i] <- sum(surv1[ft1>=ft1[i]]*diff(c(ft1[ft1>=ft1[i]], tau)))^2*nd1[i]/nrisk1[i]
    }else{
      V1i[i] <- sum(surv1[ft1>=ft1[i]]*diff(c(ft1[ft1>=ft1[i]], tau)))^2*nd1[i]/nrisk1[i]/(nrisk1[i]-nd1[i])
      }

  }


  V0 <- sum(V0i)

  V1 <- sum(V1i)



  data.frame(mu0, se0=sqrt(V0), mu1, se1=sqrt(V1), delta=delta, se_d=sqrt(V0+V1))

}
