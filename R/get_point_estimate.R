
get_point_estimate=function(trt,cox_event,cox_censor,data,M=1000,seed=NULL,cpp=TRUE,
                             memory=1024*16,local=FALSE,clmq=FALSE){

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

    sim_dt=Q(fx_clsmq_simcoun, i_bh=c(1,1,2,2), j_surv_cond=1:4,M=M,cpp=cpp, n_jobs=4,memory = memory,seed=seed,
             export = list(bh_all=bh_all,
                           s_all=s_all,
                           cpp=cpp,
                           M=M,
                           simulate_counterfactuals=simulate_counterfactuals),pkgs = c('survival','Rcpp'))

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

