stan_mod <- stan_model("inst/stan/seir.stan")

cases <- c(
  0, 0, 1, 3, 1, 8, 0, 6, 5, 0, 7, 7, 18, 9, 22, 38, 53, 45, 40,
  77, 76, 48, 67, 78, 42, 66, 67, 92, 16, 70, 43, 53, 55, 53, 29,
  26, 37, 25, 45, 34, 40, 35
)

s1 <- create_segments_vector("2020-01-01", length(cases),
  segments = c("2020-01-13"), values = c(0.1, 0.2)
)

transmission_vec <- create_ramp_vector("2020-01-01", length(cases),
  start_ramp = "2020-02-12",
  ramp_length = 7 * 6, ramp_max = 1.5
)

voc_pars <- list(
  "present" = TRUE,
  "transmission_vec" = transmission_vec
)

m <- fit_seir(
  cases,
  chains = 1,
  stan_model = stan_mod,
  iter = 100,
  fit_type = "optimizing",
  samp_frac_fixed = s1,
  voc_pars = voc_pars
)



test_that("Wrong length of transmission_vec throws error", {
  transmission_vec <- create_ramp_vector("2020-01-01", length(cases) + 20,
    start_ramp = "2020-02-12",
    ramp_length = 7 * 4, ramp_max = 10.0
  )



  expect_error(project_seir(m,
    stan_model = stan_mod,
    forecast_days = 30,
    transmission_vec = transmission_vec
  ))
})
