get_rmst_var=function(fit,data,tau,type,n.boot=1000,parallel='clmq'){

  if(type=='boot'){

    out=boot(data=data,rmst_boot_fx,R=n.boot,fit=fit,tau=tau)

    se=sd(out$t)
    lb=quantile(out$t, prob=.025)
    ub=quantile(out$t, prob=.975)

    output=c(SE=se,lb,ub)

  }



  return(output)
}

rmst_boot_fx=function(data,idx,fit,tau){
  fit_tmp <- update(fit,data=data[idx,])
  get_rmst_est(fit_tmp,data=data[idx,],tau=tau)
}


