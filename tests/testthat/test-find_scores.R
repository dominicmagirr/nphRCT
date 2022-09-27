set.seed(1)
sim_data <- sim_events_delay(
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
  n_c=50,
  n_e=50,
  max_cal_t = 36
)

df_scores_mw<-find_scores(formula=Surv(event_time,event_status)~group,
                          data=sim_data,
                          method="mw",
                          t_star = 4
)
df_scores_rmst<-find_scores(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="rmst",
                            tau = 4
)
df_scores_ms<-find_scores(formula=Surv(event_time,event_status)~group,
                          data=sim_data,
                          method="ms",
                          tau = 4
)


save_file <- function(code){
  path <- tempfile(fileext = ".RDS")
  saveRDS(code,file = path)
  path
}

test_that("scores (mw) example snapshot", {
  expect_snapshot_file(save_file(df_scores_mw), "df_scores_mw.RDS",cran = FALSE)
})

test_that("scores (rmst) example snapshot", {
  expect_snapshot_file(save_file(df_scores_rmst), "df_scores_rmst.RDS",cran = FALSE)
})

test_that("scores (ms) example snapshot", {
  expect_snapshot_file(save_file(df_scores_ms), "df_scores_ms.RDS",cran = FALSE)
})

