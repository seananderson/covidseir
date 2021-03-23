test_that("Segments and values differ in length gives error", {
  start_date <- "2020-01-01"
  segments <- c("2020-02-01","2020-03-01")
  values <- c(1,2)
  expect_error(create_segments_vector(start_date,60,segments,values=values))
})

test_that("Check vector length and number of values match",{
  start_date <- "2020-01-01"
  segments <- c("2020-01-05","2020-01-20")
  total_days <- 30
  values <- c(0,0.5,1)
  vec <- create_segments_vector(start_date,total_days,segments,values=values)

  expect_equal(length(vec),total_days)
  expect_equal(sum(vec == 0), 4)
  expect_equal(sum(vec == 0.5), 15)
  expect_equal(sum(vec == 1), 11)
})

test_that("Check creates ramp of correct length",{
  start_date <- "2020-01-01"
  start_ramp <- "2020-01-05"
  total_days <- 30
  ramp_length <- 10
  ramp_max <- 2
  vec <- create_ramp_vector(start_date,total_days,start_ramp,ramp_length,
                         ramp_max)
  expect_equal(length(vec),total_days)
})

test_that("Vector of correct length if ramp longer than total_days",{
  start_date <- "2020-01-01"
  start_ramp <- "2020-01-05"
  total_days <- 10
  ramp_length <- 10
  ramp_max <- 2
  vec <- create_ramp_vector(start_date,total_days,start_ramp,ramp_length,
                                ramp_max)
  expect_equal(length(vec),total_days)
})
