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

save_pdf <- function(plot, width = 400, height = 400) {
  path <- tempfile(fileext = ".pdf")
  pdf(file = path, width = 400, height = 400)
  on.exit(dev.off())
  plot
  path
}

test_that("plot", {
  expect_snapshot_file(save_pdf(plot.wlrt_score(df_scores)), "plot_wlrt_score.png")
})
