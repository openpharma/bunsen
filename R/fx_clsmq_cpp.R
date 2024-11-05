fx.clsmq.cpp=function(i){

  tmp <- sample(1:nrow(data), size=nrow(data), replace=TRUE)

  tmp.dt <- data[tmp, ]

  cox_event_tmp <- update(cox_event,data=tmp.dt)

  cox_censor_tmp <- update(cox_censor,data=tmp.dt)

  hr_tmp <- get_point_estimate(trt=trt,cox_event = cox_event_tmp,cox_censor = cox_censor_tmp,data=tmp.dt,M=M,seed=seed,cpp=cpp,
                               memory=memory,local=local,clmq=clmq)
  return(hr_tmp)
}

