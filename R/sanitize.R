#' Check the coxph model fit and model specifications.
#'
#' Check the coxph model fit including model class and covariates.
#' See if the model is supported for estimating the marginal treatment effects.
#'
#' @param model A coxph model from survival package.
#' @param ... Parameters for other methods.
#' @keywords internal
#' @export
sanitize_coxmodel <- function(model, ...) {
  UseMethod("sanitize_coxmodel")
}

#' Check if the model is supported.
#'
#' At the moment, only coxph is supported for the time-to-event endpoints.
#'
#' @param model A coxph model from survival package.
#' @param ... Parameters for other methods.
#' @keywords internal
#' @export
sanitize_coxmodel.default <- function(model, ...) {
  if (!inherits(model, "coxph")) {
    msg <- c(sprintf(
      'model of class "%s" is not supported and please use coxph instead.',
      class(model)[1]
    ))
    stop(msg, call. = FALSE)
  }
}


#' Check if the coxph model is correctly specified.
#'
#' Check the covariates of the coxph model.
#'
#' @param model A coxph model from survival package.
#' @param trt Character. Name of the treatment assignment variable.
#' @keywords internal
#' @export
sanitize_coxmodel.coxph <- function(model, trt, ...) {
  # check strata

  if (any(grepl("strata", attr(model$terms, "term.labels")))) stop(sprintf(c("strata is not supported at the moment.")), call. = FALSE)

  # check if trt included in the model

  if (!trt %in% attr(model$terms, "term.labels")) stop(sprintf(c("treatment variable ", trt, " is not included in the model.")), call. = FALSE)

  # check trt-covariate interactions

  if (isTRUE(grepl(trt, attr(model$terms, "term.labels")[attr(model$terms, "order") > 1]))) stop("treatment-covariate interaction is not supported at the moment.", call. = FALSE)

  # check trt levels
  data <- model.frame(model)

  if (length(na.omit(unique(data[, trt]))) != 2) stop(sprintf(c("treatment variable ", trt, " must have 2 levels.")), call. = FALSE)

}
