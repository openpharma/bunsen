get_marginal_effect=function(trt, cox_event, cox_censor, data, M, SE=TRUE, seed=NULL, cpp=TRUE,
                             memory=1024*32,local_se=FALSE,clmq_se=FALSE,clmq_hr=TRUE,clmq_local=TRUE,
                             n.boot=1000,n_jobs=100,local_cores=1,local_memory=1024*32){

  if(length(na.omit(unique(data[,trt])))!=2) stop(sprintf(c("treatment variable ",trt," must have 2 levels.")), call. = FALSE)

  sanitize_coxmodel(cox_event,trt)
  sanitize_coxmodel(cox_censor,trt)


  output=get_point_estimate(trt,cox_event,cox_censor,data,M,seed,cpp,
                        memory=local_memory,local=clmq_local,clmq=clmq_hr)

  if(SE){
    SE=get_variance_estimation(trt,data,M,n.boot,n_jobs,memory,cpp,local=local_se,clmq=clmq_se,seed,local_cores,local_memory)
    output=c(beta=output,SE)
  }

  return(output)
}
