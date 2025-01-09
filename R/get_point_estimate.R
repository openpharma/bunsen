#'Calculate the marginal treatment effects (only the point estimate) for hazard ratio (HR) adjusting covariates in clinical trials
#'
#'@description
#'This function only estimates the marginal logHR. For a complete estimation including the standard error, see \link[bunsen]{get_marginal_effect}.
#'
#'@details
#'Use the simulation approach from Daniel et al. (2020) to estimate the marginal HR when adjusting covariates.
#'This function uses the parallel computation via remote LSF and local multiprocess. Additional
#'feature also includes the C++ optimization that can speed up the calculation.
#'
#'@param trt Character. Variable name of the treatment assignment. Only support two arm trial at the moment.
#'@param cox_event Object. A coxph model using the survival time and survival status.
#'@param cox_censor Object. A coxph model using the survival time and 1-survival status.
#'@param data A data frame used for cox_event and cox_censor.
#'@param M Numeric. The number of simulated counterfactual patients. Suggest to set above 1,000,000 to get robust estimation but it is time comsuming,
#'@param seed Numeric. Random seed for simulation.
#'@param cpp Bool. True for using C++ optimization. False for not using C++ optimization. This requires cpp package installed.
#'@param memory Numeric. Memory allocation for the remote workers.
#'@param clmq Bool. True for calculating point estimate (marginal HR) using parallel computation via clustermq.
#'@param local Bool. True for calculating point estimate (marginal HR) using local multiprocess in remote workers. This is only useful when clmq_hr = TRUE.
#'@param local_cores Numeric. Number of cores or processes used in local multiprocess. This is only useful when local = TRUE. Default = 1.
#'
#'@return The marginal beta (logHR)
#'
#'@references
#'Daniel R, Zhang J, Farewell D. Making apples from oranges:
#'Comparing noncollapsible effect estimators and their standard errors
#'after adjustment for different covariate sets.
#'Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#'@export
#'@importFrom survival basehaz
#'@importFrom clustermq Q
#'@importFrom Rcpp sourceCpp
#'@examples
#'library(survival)
#'data('oak')
#'
#'cox_event <- coxph(Surv(OS, os.status) ~ trt+btmb+pdl1, data=oak)
#
#'cox_censor <- coxph(Surv(OS, 1-os.status) ~trt+btmb+pdl1, data=oak)
#
#'get_point_estimate(trt = 'trt',cox_event=cox_event,cox_censor,M=1000,data=oak,seed = 1)




get_point_estimate=function(trt,cox_event,cox_censor,data,M=1000,seed=NULL,cpp=TRUE,
                             memory=1024*16,local=FALSE,clmq=FALSE,local_cores=1){

  bh <- basehaz(cox_event, centered = FALSE)

  bh_c <- basehaz(cox_censor, centered = FALSE)

  s1_condi=calculate_statistics(model=cox_event,data=data,type=1,trt=trt,bh=bh)
  s0_condi=calculate_statistics(model=cox_event,data=data,type=0,trt=trt,bh=bh)

  s1_condi_c=calculate_statistics(model=cox_censor,data=data,type=1,trt=trt,bh=bh_c)
  s0_condi_c=calculate_statistics(model=cox_censor,data=data,type=0,trt=trt,bh=bh_c)

  if(clmq){
    bh_all=list(bh=bh,bh_c=bh_c)
    s_all=list(s1_condi,s0_condi,s1_condi_c,s0_condi_c)

    if(is.null(seed)) seed=Sys.time()

    if(local) options(clustermq.scheduler="multiprocess")
    cat('Calculating point estimate in local clustermq using multiprocess...\n')

    sim_dt=Q(.fx_clsmq_simcoun, i_bh=c(1,1,2,2), j_surv_cond=1:4,M=M,cpp=cpp, n_jobs=4,memory = memory,seed=seed,
             export = list(bh_all=bh_all,
                           s_all=s_all,
                           cpp=cpp,
                           M=M,
                           simulate_counterfactuals=simulate_counterfactuals),template = list(cores = local_cores),pkgs = c('survival','Rcpp'))

    output=calculate_trt_effect(sim_out_1d =sim_dt[[1]] ,sim_out_0d = sim_dt[[2]],sim_out_1c =sim_dt[[3]] ,sim_out_0c = sim_dt[[4]])

  }else{
    sourceCpp('cpp_functions.cpp')
    sim_out_1d=simulate_counterfactuals(bh=bh,surv_cond = s1_condi,cpp=cpp,M=M,loadcpp=FALSE)
    sim_out_0d=simulate_counterfactuals(bh=bh,surv_cond = s0_condi,cpp=cpp,M=M,loadcpp=FALSE)
    sim_out_1c=simulate_counterfactuals(bh=bh_c,surv_cond = s1_condi_c,cpp=cpp,M=M,loadcpp=FALSE)
    sim_out_0c=simulate_counterfactuals(bh=bh_c,surv_cond = s0_condi_c,cpp=cpp,M=M,loadcpp=FALSE)
    output=calculate_trt_effect(sim_out_1d,sim_out_0d,sim_out_1c,sim_out_0c)
  }

  return(output)

}


.fx_clsmq_simcoun=function(i_bh,j_surv_cond,M,cpp){
  out=simulate_counterfactuals(bh=bh_all[[i_bh]],surv_cond = s_all[[j_surv_cond]],cpp=cpp,M=M)
  out['bh']=i_bh
  out['surv_cond']=j_surv_cond
  return(out)
}

