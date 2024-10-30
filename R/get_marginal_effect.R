
get_marginal_effect=function(trt,cox_event,cox_censor,data,M=1000,seed=NULL,boot.SE=FALSE,n.boot=1000,n_jobs=100,cpp=TRUE){

  # if(cpp){
  #   cat('Loading C++ optimization...\n')
  #   sourceCpp('cpp_functions.cpp')
  # }

  # if(!is.null(seed)) set.seed(seed)

  bh <- basehaz(cox_event, centered = FALSE)

  bh_c <- basehaz(cox_censor, centered = FALSE)

  s1_condi=calculate_statistics(model=cox_event,data=data,type=1,trt=trt,bh=bh)
  s0_condi=calculate_statistics(model=cox_event,data=data,type=0,trt=trt,bh=bh)

  s1_condi_c=calculate_statistics(model=cox_censor,data=data,type=1,trt=trt,bh=bh_c)
  s0_condi_c=calculate_statistics(model=cox_censor,data=data,type=0,trt=trt,bh=bh_c)


  sim_out_1d=mirai({
    if(!is.null(seed)) set.seed(seed)
    require(Rcpp)
    require(survival)
    sourceCpp('cpp_functions.cpp')
    simulate_counterfactuals(bh=bh,surv_cond = s1_condi,cpp = cpp,M=M)
  },.args = list(M,bh,s1_condi,cpp,simulate_counterfactuals,seed))

  sim_out_0d=mirai({
    if(!is.null(seed)) set.seed(seed)
    require(Rcpp)
    require(survival)
    sourceCpp('cpp_functions.cpp')
    simulate_counterfactuals(bh=bh,surv_cond = s0_condi,cpp = cpp,M=M)
  },.args = list(M,bh,s0_condi,cpp,simulate_counterfactuals,seed))

  sim_out_1c=mirai({
    if(!is.null(seed)) set.seed(seed)
    require(Rcpp)
    require(survival)
    sourceCpp('cpp_functions.cpp')
    simulate_counterfactuals(bh=bh,surv_cond = s1_condi_c,cpp = cpp,M=M)
  },.args = list(M,bh,s1_condi_c,cpp,simulate_counterfactuals,seed))

  sim_out_0c=mirai({
    if(!is.null(seed)) set.seed(seed)
    require(Rcpp)
    require(survival)
    sourceCpp('cpp_functions.cpp')
    simulate_counterfactuals(bh=bh,surv_cond = s0_condi_c,cpp = cpp,M=M)
  },.args = list(M,bh,s0_condi_c,cpp,simulate_counterfactuals,seed))

  while(unresolved(sim_out_1d)|unresolved(sim_out_0d)|unresolved(sim_out_1c)|unresolved(sim_out_0c)){
    if(boot.SE){
      out_se=get_variance_estimation(fx.clsmq.cpp,cox_event,cox_censor,M,data,n.boot,n_jobs=n_jobs)
    }
  }

  out=get_point_estimate(sim_out_1d,sim_out_0d,sim_out_1c,sim_out_0c)
  output=out
  if(boot.SE) output=c(coef=output,out_se)


  # system.time({
  #   set.seed(1)
  #   sim_out_1d1=simulate_counterfactuals(bh=bh,surv_cond = s1_condi,cpp = cpp,M=M)
  #   # set.seed(2)
  #   sim_out_0d1=simulate_counterfactuals(bh=bh,surv_cond = s0_condi,cpp = cpp,M=M)
  #   # set.seed(3)
  #   sim_out_1c1=simulate_counterfactuals(bh=bh_c,surv_cond = s1_condi_c,cpp = cpp,M=M)
  #   # set.seed(4)
  #   sim_out_0c1=simulate_counterfactuals(bh=bh_c,surv_cond = s0_condi_c,cpp = cpp,M=M)
  #
  # })

  return(output)

}
