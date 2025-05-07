test_that("coxph model class is accepted", {
  cox_model <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
  expect_silent(sanitize_coxmodel(cox_model, "trt"))
})

test_that("glm model class throws error", {
  glm_model <- glm(os.status ~ trt + btmb + pdl1, data = oak, family = binomial)
  expect_error(sanitize_coxmodel(glm_model, "trt"))
})

test_that("trt variable not in model throws error", {
  cox_model <- coxph(Surv(OS, os.status) ~ trt + btmb + pdl1, data = oak)
  expect_error(sanitize_coxmodel(cox_model, "treatment"))
})

test_that("trt by covariate interactions not allowed in model", {
  cox_model <- coxph(Surv(OS, os.status) ~ trt * pdl1, data = oak)
  expect_error(sanitize_coxmodel(cox_model, "trt"))
})


# TODO
# - how to handle missing values
# - expected type of treatment, should it be a factor?
# - how many levels of treatment are supported?
# - detect whether cox model converged and warn user if not?
