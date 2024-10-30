get_variance_estimation=function(fx.clsmq.cpp,cox_event,cox_censor,M,data,n.boot,n_jobs,memory = 1024*24
){

  cat(paste0('Calculating SE in clustermq using bootstrap N = ',n.boot,'...\n'))

  out=Q(fun = fx.clsmq.cpp,i=1:n.boot,n_jobs = n_jobs,memory = memory,
        export = list(data=data,
                      trt='trt',
                      cox_event=cox_event,
                      cox_censor=cox_censor,
                      M=M,cpp=cpp,
                      get_marginal_effect=get_marginal_effect,
                      calculate_statistics=calculate_statistics,
                      get_point_estimate=get_point_estimate,
                      simulate_counterfactuals=simulate_counterfactuals,
                      boot.SE=FALSE

        ),pkgs = c('survival','Rcpp','mirai'))
  hr_se=do.call(c,out)
  se=sd(hr_se)
  lb=quantile(hr_se, prob=.025)
  ub=quantile(hr_se, prob=.975)
  output=c(se=se,lb,ub)
}


fx.clsmq.cpp=function(i){

  tmp <- sample(1:nrow(data), size=nrow(data), replace=TRUE)

  tmp.dt <- data[tmp, ]

  cox_event_tmp <- update(cox_event,data=tmp.dt)

  cox_censor_tmp <- update(cox_censor,data=tmp.dt)

  hr_tmp <- get_marginal_effect(trt,cox_event=cox_event_tmp, cox_censor = cox_censor_tmp, data=tmp.dt, M=M,boot.SE=FALSE)
  return(hr_tmp)
}
