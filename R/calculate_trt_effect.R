calculate_trt_effect=function(sim_out_1d,sim_out_0d,sim_out_1c,sim_out_0c){
    event1=data.frame(eventtime1_c=sim_out_1c$eventtime,
                      eventtime1_d=sim_out_1d$eventtime,
                      status1_d=sim_out_1d$status,
                      trt=1)
    event1$eventtime=ifelse(event1$eventtime1_d<event1$eventtime1_c,event1$eventtime1_d,event1$eventtime1_c)
    event1$status=ifelse(event1$eventtime1_d<event1$eventtime1_c,event1$status1_d,0)


    event0=data.frame(eventtime0_c=sim_out_0c$eventtime,
                      eventtime0_d=sim_out_0d$eventtime,
                      status0_d=sim_out_0d$status,
                      trt=0)
    event0$eventtime=ifelse(event0$eventtime0_d<event0$eventtime0_c,event0$eventtime0_d,event0$eventtime0_c)
    event0$status=ifelse(event0$eventtime0_d<event0$eventtime0_c,event0$status0_d,0)

    newdata <- rbind(event1[,c('eventtime','status','trt')],event0[,c('eventtime','status','trt')])

    cox_fit_new <- coxph(Surv(eventtime, status) ~ trt, newdata)

    output=as.numeric(coef(cox_fit_new))

  return(output)

}
