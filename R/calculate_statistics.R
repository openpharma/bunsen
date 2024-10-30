
calculate_statistics=function(model,data,type,trt,bh){

  dt=model.frame(model,data)[names(coefficients(model))]

  if(type==1) dt[trt]=1
  if(type==0) dt[trt]=0



  pred=predict(model,newdata = dt,type='risk',reference = 'zero')
  survf=exp(-pred%*%t(bh$hazard))

  survf_mean=colMeans(survf)
  survf_1=c(1,survf_mean[1:(length(survf_mean)-1)])
  surv_cond=survf_mean/survf_1
  return(surv_cond)
}

