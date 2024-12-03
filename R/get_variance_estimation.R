get_variance_estimation=function(trt,data,M,n.boot,n_jobs,memory=1024*32,cpp=TRUE,local=TRUE,clmq=TRUE,
                                 seed=NULL,local_cores=1,local_memory=1024*16
){
  if(is.null(seed)) seed=Sys.time()
  cat(paste0('Calculating SE in clustermq using bootstrap N = ',n.boot,'...\n'))

  options(clustermq.scheduler="LSF")


  out=Q(fun = fx.clsmq.cpp,i=1:n.boot,n_jobs = n_jobs,memory = memory,seed=seed,
        export = list(data=data,memory = local_memory,
                      cox_event=cox_event,cox_censor=cox_censor,
                      trt=trt,local=local,clmq=clmq,
                      M=M,cpp=cpp,seed=seed,
                      get_point_estimate=get_point_estimate,
                      calculate_statistics=calculate_statistics,
                      calculate_trt_effect=calculate_trt_effect,
                      simulate_counterfactuals=simulate_counterfactuals,
                      fx_clsmq_simcoun=fx_clsmq_simcoun
        ),template = list(cores = local_cores),pkgs = c('survival','Rcpp','clustermq'))
  hr_se=do.call(c,out)
  se=sd(hr_se)
  lb=quantile(hr_se, prob=.025)
  ub=quantile(hr_se, prob=.975)
  output=c(se=se,lb,ub)
}


