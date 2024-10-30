simulate_counterfactuals=function(bh,surv_cond,M=M,cpp=TRUE){
  if(cpp){
    newd= rbinom_matrix_vec(nrows=nrow(bh), ncols=M, surv_cond)
    index_event_d = firstZeroIndex(newd)

  }else{
    newd=matrix(rbinom(nrow(bh) * M, 1, surv_cond), nrow = nrow(bh), ncol = M)
    index_event_d = apply(newd, 2, function(x) which(x == 0)[1])
  }

  eventtime_d <- bh$time[index_event_d]

  eventtime_d[which(is.na(index_event_d))] <- max(bh$time)

  status_d <- rep(1, M)

  status_d[which(is.na(index_event_d))] <- 0
  return(list(eventtime=eventtime_d,status=status_d))
}
