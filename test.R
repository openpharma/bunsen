data('oak')

cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)

cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)

get_marginal_effect(trt = 'trt',cox_event,cox_censor,M=10000,data=oak)

use_r('get_rmst_estimate')




tau=26

fit3 <- RMST_stratCox(time=oak$OS, status=oak$os.status, arm=oak$trt, covariates=cbind(oak$btmb, oak$pdl1), tau=tau)












