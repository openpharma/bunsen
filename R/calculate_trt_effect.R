#' Calculate the marginal treatment effect using counterfactual simulations
#' @description
#' This function is used to calculate the logHR using cox model after obtaining the counterfactural simulations for potential outcomes from \link[bunsen]{simulate_counterfactuals}.
#' @details
#' The event indicator from this function uses the equation 10', I(Y0<Y0_cens) or I(Y1<Y1_cens).
#'
#' @param sim_out_1d List. A list from \link[bunsen]{simulate_counterfactuals} for cox_event in treatment group, e.g. the cox model using OS.
#' @param sim_out_0d List. A list from \link[bunsen]{simulate_counterfactuals} for cox_event in control group, e.g. the cox model using OS.
#' @param sim_out_1c List. A list from \link[bunsen]{simulate_counterfactuals} for cox_event in treatment group, e.g. the cox model using 1-OS.
#' @param sim_out_0c List. A list from \link[bunsen]{simulate_counterfactuals} for cox_event in control group, e.g. the cox model using 1-OS.
#'
#' @return The marginal beta (logHR)
#'
#' @references
#' Daniel R, Zhang J, Farewell D. Making apples from oranges:
#' Comparing noncollapsible effect estimators and their standard errors
#' after adjustment for different covariate sets.
#' Biom J. 2021;63(3):528-557. doi:10.1002/bimj.201900297
#' @export
#' @importFrom survival coxph
#' @examples
#' library(survival)
#' data("oak")
#'
#' cox_event <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
#' #
#' cox_censor <- coxph(Surv(OS, 1 - os.status) ~ trt + btmb + pdl1, data = oak)
#' bh <- basehaz(cox_event, centered = FALSE)
#' bh_c <- basehaz(cox_censor, centered = FALSE)
#' s_condi <- calculate_statistics(model = cox_event, trt = "trt")
#' s_condi_c <- calculate_statistics(model = cox_censor, trt = "trt")
#' sim_out_1d <- simulate_counterfactuals(bh = bh, surv_cond = s_condi$surv_cond1, cpp = FALSE, M = 1000)
#' sim_out_0d <- simulate_counterfactuals(bh = bh, surv_cond = s_condi$surv_cond0, cpp = FALSE, M = 1000)
#' sim_out_1c <- simulate_counterfactuals(bh = bh_c, surv_cond = s_condi_c$surv_cond1, cpp = FALSE, M = 1000)
#' sim_out_0c <- simulate_counterfactuals(bh = bh_c, surv_cond = s_condi_c$surv_cond0, cpp = FALSE, M = 1000)
#'
#' output <- calculate_trt_effect(sim_out_1d, sim_out_0d, sim_out_1c, sim_out_0c)
#'

calculate_trt_effect <- function(
  sim_out_1d,
  sim_out_0d,
  sim_out_1c,
  sim_out_0c
) {
  event1 <- data.frame(
    eventtime1_c = sim_out_1c$eventtime,
    eventtime1_d = sim_out_1d$eventtime,
    status1_d = sim_out_1d$status,
    trt = 1
  )
  event1$eventtime <- ifelse(
    event1$eventtime1_d < event1$eventtime1_c,
    event1$eventtime1_d,
    event1$eventtime1_c
  )
  event1$status <- ifelse(
    event1$eventtime1_d < event1$eventtime1_c,
    event1$status1_d,
    0
  )

  event0 <- data.frame(
    eventtime0_c = sim_out_0c$eventtime,
    eventtime0_d = sim_out_0d$eventtime,
    status0_d = sim_out_0d$status,
    trt = 0
  )

  event0$eventtime <- ifelse(
    event0$eventtime0_d < event0$eventtime0_c,
    event0$eventtime0_d,
    event0$eventtime0_c
  )
  event0$status <- ifelse(
    event0$eventtime0_d < event0$eventtime0_c,
    event0$status0_d,
    0
  )

  newdata <- rbind(
    event1[, c('eventtime', 'status', 'trt')],
    event0[, c('eventtime', 'status', 'trt')]
  )

  cox_fit_new <- coxph(Surv(eventtime, status) ~ trt, newdata)


  output <- as.numeric(coef(cox_fit_new))

  return(output)
}
