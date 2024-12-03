sanitize_coxmodel = function(object, ...) {
  UseMethod("sanitize_coxmodel")
}



sanitize_coxmodel.default <- function(model,...) {
  if (!inherits(model, "coxph")) {
    msg <- c(sprintf('model of class "%s" is not supported and please use coxph instead.', class(model)[1]))
    stop(msg, call. = FALSE)
  }

}


sanitize_coxmodel.coxph=function(model,trt){

  # check strata

  if(any(grepl('strata', attr(model$terms,'term.labels')))) stop(sprintf(c("strata is not supported at the moment.")), call. = FALSE)

  # check if trt included in the model

  if(!trt%in%attr(model$terms,'term.labels')) stop(sprintf(c("treatment variable ",trt," is not included in the model.")), call. = FALSE)

  # check trt-covariate interactions

  if(grepl(trt, attr(model$terms,'term.labels')[attr(model$terms,'order')>1])) stop("treatment-covariate interaction is not supported at the moment.", call. = FALSE)




}

