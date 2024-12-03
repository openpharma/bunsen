library(survival)
library(boot)
library(clustermq)
data('oak')
# oak[5,c(2)]=NA
# oak[10,c(3)]=NA
# oak[6,c(1)]=NA


cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)

cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)

get_marginal_effect(trt = 'trt',cox_event=cox_event,cox_censor,SE=FALSE,M=10000,data=oak,seed = 1)


RMST_stratCox(time=oak$OS, status=oak$os.status, arm=oak$trt, covariates=cbind(oak$btmb, oak$pdl1), tau=tau)


data('oak')

tau=26

fit_new <- coxph(Surv(OS, os.status) ~ btmb+strata(trt)+pdl1, data=oak)

get_rmst_est(fit=fit,data=oak,tau=tau)

set.seed(960159)
a=get_rmst_var(fit=fit,data=oak,tau=tau,type = 'boot',n.boot = 1000)


tau=26
time=oak$OS
status=oak$os.status
arm=oak$trt
covariates=cbind(oak$btmb, oak$pdl1)

get_rmst_est(time, status, arm, covariates,tau,SE='delta')


get_rmst_est(time, status, arm, covariates,tau,SE='boot',n.boot = 1000)
