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

save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  png(path, width = width, height = height)
  on.exit(dev.off())
  code
  
  path
}

expect_snapshot_file(save_png(plot(df_scores_mw)), paste0("plot_mw.png"),cran = TRUE)
expect_snapshot_file(save_png(plot(df_scores_rmst)), paste0("plot_rmst.png"),cran = TRUE)
expect_snapshot_file(save_png(plot(df_scores_ms)), paste0("plot_ms.png"),cran = TRUE)
