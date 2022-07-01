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


df<-df_scores$df
test_that("should add to zero by definition", {
  expect_equal(sum(c(rowSums(df[,grep("n_event_",names(df))])*df$score_event+
                 rowSums(df[,grep("n_censored_",names(df))])*df$score_censored)), 0)
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
