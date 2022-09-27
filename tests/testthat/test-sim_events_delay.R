set.seed(1)

#power model example
sim_data_power <- sim_events_delay(
  event_model=list(
    duration_c = 36,
    duration_e = c(6,30),
    lambda_c = log(2)/9,
    lambda_e = c(log(2)/9,log(2)/18)
  ),
  recruitment_model=list(
    rec_model="power",
    rec_period = 12,
    rec_power = 1
  ),
  n_c=5,
  n_e=5,
  max_cal_t = 36
)

test_that("returns a data frame", {expect_s3_class(sim_data_power, "data.frame")})
test_that("correct number rows", {expect_identical(nrow(sim_data_power), 10L)})
test_that("correct number columns", {expect_identical(ncol(sim_data_power), 3L)})
test_that("no NAs", {expect_false(any(is.na(sim_data_power)))})

save_file <- function(code){
  path <- tempfile(fileext = ".RDS")
  saveRDS(code,file = path)
  path
}

test_that("power model example snapshot", {
  expect_snapshot_file(save_file(sim_data_power), "sim_data_power.RDS",cran = FALSE)
})


#piecewise constant model
sim_data_pw <- sim_events_delay(
  event_model=list(
    duration_c = 36,
    duration_e = c(6,30),
    lambda_c = log(2)/9,
    lambda_e = c(log(2)/9,log(2)/18)
  ),
  recruitment_model=list(
    rec_model="pw_constant",
    rec_rate = seq(1,2,length=6),
    rec_duration = rep(6,6)
  ),
  n_c=5,
  n_e=5,
  max_cal_t = 36
)

test_that("returns a data frame", {expect_s3_class(sim_data_pw, "data.frame")})
test_that("correct number rows", {expect_identical(nrow(sim_data_pw), 10L)})
test_that("correct number columns", {expect_identical(ncol(sim_data_pw), 3L)})
test_that("no NAs", {expect_false(any(is.na(sim_data_pw)))})

test_that("piecewise constant model example snapshot", {
  expect_snapshot_file(save_file(sim_data_pw), "sim_data_pw.RDS",cran = FALSE)
})

#Constant
sim_data_constant <- sim_events_delay(
  event_model=list(
    duration_c = 36,
    duration_e = c(6,30),
    lambda_c = log(2)/9,
    lambda_e = c(log(2)/9,log(2)/18)
  ),
  recruitment_model=list(
    rec_model="pw_constant",
    rec_rate = 2,
    rec_duration = 36
  ),
  n_c=5,
  n_e=5,
  max_cal_t = 36
)

test_that("returns a data frame", {expect_s3_class(sim_data_constant, "data.frame")})
test_that("correct number rows", {expect_identical(nrow(sim_data_constant), 10L)})
test_that("correct number columns", {expect_identical(ncol(sim_data_constant), 3L)})
test_that("no NAs", {expect_false(any(is.na(sim_data_constant)))})

test_that("piecewise constant model example snapshot", {
  expect_snapshot_file(save_file(sim_data_constant), "sim_data_constant.RDS",cran = FALSE)
})


#Errors
test_that("rec_rate is negative error", {
  expect_error(sim_events_delay(
    event_model=list(
      duration_c = 36,
      duration_e = c(6,30),
      lambda_c = log(2)/9,
      lambda_e = c(log(2)/9,log(2)/18)
    ),
    recruitment_model=list(
      rec_model="pw_constant",
      rec_rate = seq(-1,2,length=6),
      rec_duration = rep(6,6)
    ),
    n_c=5,
    n_e=5,
    max_cal_t = 36
  ),
  "rec_rate should be non-negative")
})

test_that("durations are correct", {
  expect_error(sim_data <- sim_events_delay(
    event_model=list(
      duration_c = 36,
      duration_e = c(6,29),
      lambda_c = log(2)/9,
      lambda_e = c(log(2)/9,log(2)/18)
    ),
    recruitment_model=list(
      rec_model="power",
      rec_period = 12,
      rec_power = 1
    ),
    n_c=5,
    n_e=5,
    max_cal_t = 36),"The total length of durations in duration_c and duration_e must not be less than max_cal_t")})
