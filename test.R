
# package setup -----------------------------------------------------------
library(usethis)
use_news_md()
use_rcpp()
use_package('survival')
usethis::use_package("clustermq")
usethis::use_package("Rcpp")
usethis::use_package_doc()
usethis::use_roxygen_md()

# test code ---------------------------------------------------------------


library(survival)
library(clustermq)
data('oak')

cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)

cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)

get_marginal_effect(trt = 'trt',cox_event=cox_event,cox_censor,SE=TRUE,M=1000,data=oak,seed = 1)


data('oak')

tau=26

time=oak$OS

status=oak$os.status

arm=oak$trt

covariates=oak[,c('btmb','pdl1')]

get_rmst_est(time, status, arm, covariates,tau,SE='delta')

set.seed(2024)

get_rmst_est(time, status, arm, covariates,tau,SE='boot',n.boot = 1000)

library(survival)
data('oak')

cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)
cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)
bh <- basehaz(cox_event, centered = FALSE)
bh_c <- basehaz(cox_censor, centered = FALSE)
s1_condi=calculate_statistics(model=cox_event,data=oak,type=1,trt='trt',bh=bh)
s0_condi=calculate_statistics(model=cox_event,data=oak,type=0,trt='trt',bh=bh)
s1_condi_c=calculate_statistics(model=cox_censor,data=oak,type=1,trt='trt',bh=bh_c)
s0_condi_c=calculate_statistics(model=cox_censor,data=oak,type=0,trt='trt',bh=bh_c)
sim_out_1d=simulate_counterfactuals(bh=bh,surv_cond = s1_condi,cpp=FALSE,M=1000)
sim_out_0d=simulate_counterfactuals(bh=bh,surv_cond = s0_condi,cpp=FALSE,M=1000)
sim_out_1c=simulate_counterfactuals(bh=bh_c,surv_cond = s1_condi_c,cpp=FALSE,M=1000)
sim_out_0c=simulate_counterfactuals(bh=bh_c,surv_cond = s0_condi_c,cpp=FALSE,M=1000)

output=calculate_trt_effect(sim_out_1d,sim_out_0d,sim_out_1c,sim_out_0c)
library(survival)
data('oak')

tau=26
time=oak$OS
status=oak$os.status
arm=oak$trt
covariates=oak[,c("btmb","pdl1")]

dt=as.data.frame(cbind(time,status, arm, covariates))
fit=coxph(Surv(time, status) ~ btmb+pdl1+strata(arm),data = dt)
delta=rmst_point_estimate(fit,dt=dt,tau)
delta$output
