sim_data <- sim_events_delay(
  n_c = 5,
  n_e = 5,
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_period = 12,
  rec_power = 1,
  max_cal_t = 36
)
test_that("returns a data frame", {expect_s3_class(sim_data, "data.frame")})
test_that("correct number rows", {expect_identical(nrow(sim_data), 10L)})
test_that("correct number columns", {expect_identical(ncol(sim_data), 3L)})
test_that("no NAs", {expect_false(any(is.na(sim_data)))})