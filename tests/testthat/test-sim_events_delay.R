set.seed(1)
rec_c <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=5)
rec_e <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=5)
sim_data <- sim_events_delay(
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_times_c = rec_c,
  rec_times_e = rec_e,
  max_cal_t = 36
)
test_that("returns a data frame", {expect_s3_class(sim_data, "data.frame")})
test_that("correct number rows", {expect_identical(nrow(sim_data), 10L)})
test_that("correct number columns", {expect_identical(ncol(sim_data), 3L)})
test_that("no NAs", {expect_false(any(is.na(sim_data)))})


test_that("rec_times between 0 and max_cal_t", {
  expect_error(sim_events_delay(
    delay_e = 6,
    lambda_c = log(2)/9,
    lambda_e_1 = log(2)/9,
    lambda_e_2 = log(2)/18,
    rec_times_c = rec_c,
    rec_times_e = rec_e,
    max_cal_t = 11
  ),
  "patients must enroll between 0 and max_cal_t")
})
