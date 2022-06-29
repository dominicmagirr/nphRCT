sim_data <- sim_events_delay(
  n_c = 1000,
  n_e = 1000,
  delay_e = 6,
  lambda_c = log(2)/9,
  lambda_e_1 = log(2)/9,
  lambda_e_2 = log(2)/18,
  rec_period = 12,
  rec_power = 1,
  max_cal_t = 36
)
df_weights_mw<-find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="mw",
                            t_star = 4
)
df_weights_mw_s_star<-find_weights(formula=Surv(event_time,event_status)~group,
                                   data=sim_data,
                                   wlr="mw",
                                   s_star = 0.5
)
df_weights_mw_0<-find_weights(formula=Surv(event_time,event_status)~group,
                              data=sim_data,
                              wlr="mw",
                              t_star = 0
)
df_weights_cens<-find_weights(formula=Surv(event_time,event_status)~group,
                              data=sim_data,
                              wlr="mw",
                              t_star = 0,
                              include_cens=FALSE
)
df_weights_lr<-find_weights(formula=Surv(event_time,event_status)~group,
                         data=sim_data,
                         wlr="lr"
)
df_weights_fh<-find_weights(formula=Surv(event_time,event_status)~group,
                         data=sim_data,
                         wlr="fh",
                         rho = 0,
                         gamma=1
)

for (df_weights in 
     list(df_weights_mw,df_weights_lr,df_weights_fh,
          df_weights_mw_s_star,df_weights_mw_0,df_weights_cens)){
  test_that("returns a data frame", {expect_s3_class(df_weights, "data.frame")})
  
  test_that("correct number of rows", {
    expect_true(nrow(df_weights)<=2000)
  })
  
  test_that("check time is strictly increasing", {
    expect_true(all(diff(df_weights$t_j)>0))
  })
  
  
  test_that("positive number of events", {
    expect_true(all(df_weights$n_event>=0))
  })
  
  test_that("positive number of censored", {
    expect_true(all(df_weights$S_hat<=df_weights$S_hat_minus))
  })
  
  test_that("weights are increasing", {
    expect_true(all(diff(df_weights$w)>=0))
  })
  
  test_that("number at risk are decreasing", {
    expect_true(all(diff(df_weights$n_risk)<0))
  })
}

#Test errors

test_that("specifying rho with lr test", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="lr",
                            rho = 0),
               "do not specify rho or gamma for log rank test")
})

test_that("specifying t_star with lr test", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="lr",
                            t_star=4),
               "do not specify t_star or s_star for log rank test ")
})

test_that("specifying rho with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="mw",
                            rho = 0,
                            t_star=4),
               "do not specify rho or gamma for modestly weighted log rank test")
})

test_that("null gamma with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="fh",
                            rho = 0),
               "specify rho and gamma")
})

test_that("specify t_star with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="fh",
                            rho = 0,
                            gamma=1,
                            t_star=4),
               "do not specify t_star or s_star for Fleming-Harrington test")
})
test_that("negative gamma with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="fh",
                            rho = 0,
                            gamma=-1),
               "rho and gamma must be non-negative")
})

test_that("s_star not a probability", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="mw",
                            s_star = 0.5,
                            t_star= 4),
               "must specify either t_star or s_star")
})

test_that("s_star and t_star not a specified", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="mw"),
               "must specify either t_star or s_star")
})

test_that("s_star not a probability", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            wlr="mw",
                            s_star = -0.5),
               "s_star must a probability")
})

sim_data_strata<-sim_data
sim_data_strata$ecog<-rep(c(0,1),1000)
test_that("don't allow strata",
  {expect_error(find_weights(formula=Surv(event_time,event_status)~group+strata(ecog),
                            data=sim_data_strata,
                            wlr="mw",
                            s_star = 0.5), "does not account for strata")}
)
