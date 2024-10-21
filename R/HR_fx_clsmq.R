HR_fx_clsmq=function(i){

  tmp <- sample(1:nrow(data), size=nrow(data), replace=TRUE)

  tmp.dt <- data[tmp, ]

  cox_event_tmp <- update(cox_event,data=tmp.dt)

  cox_censor_tmp <- update(cox_censor,data=tmp.dt)

  hr_tmp <- marginal_HR(trt,cox_event=cox_event_tmp, cox_censor = cox_censor_tmp, data=tmp.dt, M=M)
  return(hr_tmp)
}
