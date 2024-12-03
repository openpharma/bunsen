get_rmst_var=function(fit,data,tau,type,n.boot=1000){

  if(type=='boot'){

    out=boot(data=data,rmst_boot_fx,R=n.boot,fit=fit,tau=tau)

    se=sd(out$t,na.rm = T)
    lb=quantile(out$t, prob=.025,na.rm=T)
    ub=quantile(out$t, prob=.975,na.rm=T)

    output=list(c(SE=se,lb,ub),t=out$t)

  }

  if(type=='delta'){

    dt=model.frame(fit_new,data)
    names(dt)=c('time','status',gsub("strata\\(([^\\)]+)\\)", "\\1", attr(fit_new$terms,'term.labels')))
    trt=gsub("strata\\(([^\\)]+)\\)", "\\1", attr(fit_new$terms,'term.labels'))[grepl("strata",attr(fit_new$terms,'term.labels'))]
    dt[trt]=


    dt0=
    dt1=data[data[trt]==unique(na.omit(data[,trt]))[2],]

    sigma_new=fit_new$var*fit_new$n

    xb0=na.omit(predict(fit_new,newdata = dt0,type='lp',reference = 'zero'))



  }



  return(output)
}




rmst_boot_fx=function(data,idx,fit,tau){

  fit_tmp <- update(fit,data=data[idx,])
  tmax=basehaz(fit)
  tmax0=max(tmax$time[tmax$strata==unique(tmax$strata)[1]])
  tmax1=max(tmax$time[tmax$strata==unique(tmax$strata)[2]])
  if(tmax0<tau |tmax1<tau) {
    output=NA
  }else{
    output=get_rmst_est(fit_tmp,data=data[idx,],tau=tau)
  }
  return(output)
}


