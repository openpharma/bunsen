fx_clsmq_simcoun=function(i_bh,j_surv_cond,M,cpp){
  out=simulate_counterfactuals(bh=bh_all[[i_bh]],surv_cond = s_all[[j_surv_cond]],cpp=cpp,M=M)
  out['bh']=i_bh
  out['surv_cond']=j_surv_cond
  return(out)
}
