set.seed(1)
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
df_scores<-find_scores(formula=Surv(event_time,event_status)~group,
  data=sim_data,
  wlr="mw",
  t_star = 4
)

test_that("should add to zero by definition", {
  expect_equal(sum(df_scores$df$n_event*df_scores$df$score_event)+
                 sum(df_scores$df$n_censor*df_scores$df$score_cens), 0)
})

save_png <- function(code) {
  path <- tempfile(fileext = ".png")
  png(file = path)
  on.exit(dev.off())
  code
  path
}

test_that("plot", {
  expect_snapshot_file(save_png(	
    plot.wlrt_score(df_scores)), "plot_wlrt_score.png",cran = TRUE)
})

save_table <- function(code) {
  path <- tempfile(fileext = ".png")
  png(file = path)
  on.exit(dev.off())
  code
  path
}
