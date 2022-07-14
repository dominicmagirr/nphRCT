rec_c <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=50)
rec_e <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=50)
sim_data <- sim_events_delay(
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_times_c = rec_c,
  rec_times_e = rec_e,
  max_cal_t = 36
)
sim_data_0 <- sim_data
sim_data_0$ecog=0
rec_c <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=50)
rec_e <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=50)
sim_data_1 <- sim_events_delay(
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_times_c = rec_c,
  rec_times_e = rec_e,
  max_cal_t = 36
)
sim_data_1$ecog=1
sim_data_strata<-rbind(sim_data_0,sim_data_1)
sim_data_strata_2<-cbind(sim_data_strata,sex=rep(c("M","F"),times=100))
rec_c <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=2)
rec_e <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=2)
sim_data_2 <- sim_events_delay(
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_times_c = rec_c,
  rec_times_e = rec_e,
  max_cal_t = 36
)
sim_data_strata_2<-cbind(sim_data_strata,sex=rep(c("M","F"),times=100))
sim_data_strata_small<-rbind(sim_data_0,cbind(sim_data_2,ecog=1))


save_file <- function(code) {
  path <- tempfile(fileext = ".RDS")
  saveRDS(code,file = path)
  path
}

test_that("example setting t_star", {
  expect_snapshot_file(save_file(	
    wlrt(formula=Surv(event_time,event_status)~group,
         data=sim_data,
         wlr="mw",
         t_star = 4
    )), "out_t_star.RDS",cran = TRUE)
})
test_that("example setting s_star", {
  expect_snapshot_file(save_file(	
    wlrt(formula=Surv(event_time,event_status)~group,
         data=sim_data,
         wlr="mw",
         s_star = 0.5
    )), "out_s_star.RDS",cran = TRUE)
})
test_that("example with 1 strata", {
  expect_snapshot_file(save_file(	
    wlrt(formula=Surv(event_time,event_status)~group+strata(ecog),
         data=sim_data_strata,
         wlr="mw",
         t_star = 4
    )), "out_strata.RDS",cran = TRUE)
})

test_that("example with 2 strata", {
  expect_snapshot_file(save_file(	
    wlrt(formula=Surv(event_time,event_status)~group+strata(ecog)+strata(sex),
         data=sim_data_strata_2,
         wlr="mw",
         t_star = 4
    )), "out_strata_2.RDS",cran = TRUE)
})


#Test errors

test_that("only one treatment arm indicator", {
  expect_error(wlrt(formula=Surv(event_time,event_status)~group+sex+strata(ecog),
                    data=sim_data_strata_2,
                    wlr="mw",
                    t_star = 4
  ),
  "Formula must contain only one treatment arm indicator")
})

test_that("Minimum stratum size is 5", {
  expect_error(wlrt(formula=Surv(event_time,event_status)~group+strata(ecog),
                    data=sim_data_strata_small,
                    wlr="mw",
                    t_star = 4
  ),
  "Minimum stratum size is 5")
})

