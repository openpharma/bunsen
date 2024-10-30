'================================================================'
# Functions for covairate adjustment -- UK VAN project
#
# Author: Xinlei Deng
# R Version: R version 4.3.1 (2023-06-16)
# Last Modified Date (Year-MM-DD): 2024-10-10
'================================================================'

# trt: a character value indicating the variable name of the treatment
# cox_event: cox model with event as the endpoint
# cox_censor: cox model with 1-event as the endpoint
# data: data used to fit the cox models
# M: simulated number of patients
# seed: random seed for simulation to get the same resultr
# boot.SE: using bootstrap to calculate SE
# n.boot: number of bootstrap
# n_jobs: The number of LSF jobs to submit in clustermq
#
# require(dplyr)
# require(survival)
# require(clustermq)
#
# marginal_HR=function(trt,cox_event,cox_censor,data,M=1000,seed=NULL,boot.SE=FALSE,n.boot=1000,n_jobs=100){
#
#   if(!is.null(seed)) set.seed(seed)
#
#   dt_1=dt_0=model.frame(cox_event,data)[names(coefficients(cox_event))]
#   dt_1[trt]=1
#   dt_0[trt]=0
#
#   # calculate Z1 and Z0 from event cox
#   bh <- basehaz(cox_event, centered = FALSE)
#
#   pred1=predict(cox_event,newdata = dt_1,type='risk',reference = 'zero')
#   pred0=predict(cox_event,newdata = dt_0,type='risk',reference = 'zero')
#
#
#   s1=exp(-pred1%*%t(bh$hazard))
#   s0=exp(-pred0%*%t(bh$hazard))
#
#   s1_mean=colMeans(s1)
#   s1_mean_copy=c(1,s1_mean[1:(length(s1_mean)-1)])
#   s1_condi=s1_mean/s1_mean_copy
#
#   s0_mean=colMeans(s0)
#   s0_mean_copy=c(1,s0_mean[1:(length(s0_mean)-1)])
#   s0_condi=s0_mean/s0_mean_copy
#
#   # calculate Z1 and Z0 from censor cox
#   bh_c <- basehaz(cox_censor, centered = FALSE)
#
#   pred1_c=predict(cox_censor,newdata = dt_1,type='risk',reference = 'zero')
#   pred0_c=predict(cox_censor,newdata = dt_0,type='risk',reference = 'zero')
#
#
#   s1_c=exp(-pred1_c%*%t(bh_c$hazard))
#   s0_c=exp(-pred0_c%*%t(bh_c$hazard))
#
#
#   s1_mean_c=colMeans(s1_c)
#   s1_mean_c_copy=c(1,s1_mean_c[1:(length(s1_mean_c)-1)])
#   s1_condi_c=s1_mean_c/s1_mean_c_copy
#
#   s0_mean_c=colMeans(s0_c)
#   s0_mean_c_copy=c(1,s0_mean_c[1:(length(s0_mean_c)-1)])
#   s0_condi_c=s0_mean_c/s0_mean_c_copy
#
#   newd1=sapply(1:M, function(u){
#     rbinom(nrow(bh), 1, s1_condi)
#   })
#
#   index_event1_d <- apply(newd1, 2, function(x) which(x == 0)[1])
#
#   eventtime1_d <- bh$time[index_event1_d]
#
#   eventtime1_d[which(is.na(index_event1_d))] <- max(bh$time)
#
#   status1_d <- rep(1, M)
#
#   status1_d[which(is.na(index_event1_d))] <- 0
#
#
#   newd1_c=sapply(1:M, function(u){
#     rbinom(nrow(bh_c), 1, s1_condi_c)
#   })
#
#   index_event1_c <- apply(newd1_c, 2, function(x) which(x == 0)[1])
#
#   eventtime1_c <- bh_c$time[index_event1_c]
#
#   eventtime1_c[which(is.na(index_event1_c))] <- max(bh_c$time)
#
#   status1_c <- rep(1, M)
#
#   status1_c[which(is.na(index_event1_c))] <- 0
#
#   event1=data.frame(eventtime1_c,eventtime1_d,status1_d,trt=1)
#
#   event1=mutate(event1,eventtime=ifelse(eventtime1_d<eventtime1_c,eventtime1_d,eventtime1_c),
#                 status=ifelse(eventtime1_d<eventtime1_c,status1_d,0)
#   )
#
#
#   newd0=sapply(1:M, function(u){
#     rbinom(nrow(bh), 1, s0_condi)
#   })
#
#   index_event0_d <- apply(newd0, 2, function(x) which(x == 0)[1])
#
#   eventtime0_d <- bh$time[index_event0_d]
#
#   eventtime0_d[which(is.na(index_event0_d))] <- max(bh$time)
#
#   status0_d <- rep(1, M)
#
#   status0_d[which(is.na(index_event0_d))] <- 0
#
#
#   newd0_c=sapply(1:M, function(u){
#     rbinom(nrow(bh_c), 1, s0_condi_c)
#   })
#
#   index_event0_c <- apply(newd0_c, 2, function(x) which(x == 0)[1])
#
#   eventtime0_c <- bh_c$time[index_event0_c]
#
#   eventtime0_c[which(is.na(index_event0_c))] <- max(bh_c$time)
#
#   status0_c <- rep(1, M)
#
#   status0_c[which(is.na(index_event0_c))] <- 0
#
#   event0=data.frame(eventtime0_c,eventtime0_d,status0_d,trt=0)
#
#   event0=mutate(event0,eventtime=ifelse(eventtime0_d<eventtime0_c,eventtime0_d,eventtime0_c),
#                 status=ifelse(eventtime0_d<eventtime0_c,status0_d,0)
#   )
#
#
#   newdata <- rbind(event1[,c('eventtime','status','trt')],event0[,c('eventtime','status','trt')])
#
#   cox_fit_new <- coxph(Surv(eventtime, status) ~ trt, newdata)
#
#   output=as.numeric(coef(cox_fit_new))
#
#
#   if(boot.SE){
#     cat(paste0('Calculating SE in clustermq using bootstrap N = ',n.boot,'...\n'))
#     out=Q(fun = HR_fx_clsmq,i=1:n.boot,n_jobs = n_jobs,memory = 1024*32,export = list(data=data,
#                                                                                    trt='trt',
#                                                                                    cox_event=cox_event,
#                                                                                    cox_censor=cox_censor,
#                                                                                    M=M,
#                                                                                    marginal_HR=marginal_HR),pkgs = c('survival','dplyr'))
#     hr_se=do.call(c,out)
#     se=sd(hr_se)
#     lb=quantile(hr_se, prob=.025)
#     ub=quantile(hr_se, prob=.975)
#     output=c(coef=as.numeric(coef(cox_fit_new)),se=se,lb,ub)
#   }
#
#   return(output)
#
# }
#
#


library(openxlsx)
library(dplyr)
library(survival)
library(clustermq)
library(Rcpp)
library(mirai)
# setDTthreads(4)

oakdata=read.xlsx('OAK data.xlsx',sheet = 'OAK_Clinical_Data')

oak_simdt=oakdata %>% filter(BEP=='Y',EGFRMUT!='POSITIVE',EML4MUT!='POSITIVE') %>%
  mutate(pfs.status=1-PFS.CNSR,os.status=1-OS.CNSR, orr=as.numeric(BCOR %in% c("PR", "CR")),
         btmb=as.numeric(btmb),pdl1=ifelse(TC1IC1=="UNKNOWN", NA, TC1IC1)) %>%
  mutate(pdl1=ifelse(pdl1=='TC1/2/3 or IC1/2/3',1,0),trt=ifelse(TRT01P=='MPDL3280A',1,0))

oak_simdt=oak_simdt %>% select(trt,btmb,pdl1,OS,os.status)
oak_simdt=na.omit(oak_simdt)



cox1 <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak_simdt)

cox2 <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak_simdt)

M=1000000

trt='trt'

cox_event=cox1

cox_censor=cox2
system.time({
  get_marginal_effect(trt,cox_event,cox_censor,data=oak_simdt,M=1000000,seed=NULL,boot.SE=FALSE,cpp=TRUE)
})
















