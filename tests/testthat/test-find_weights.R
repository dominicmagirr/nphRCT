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
  n_c=1000,
  n_e=1000,
  max_cal_t = 36
)
weights_mw<-find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="mw",
                            t_star = 4
)
weights_mw_s_star<-find_weights(formula=Surv(event_time,event_status)~group,
                                   data=sim_data,
                                   method="mw",
                                   s_star = 0.5
)
weights_mw_0<-find_weights(formula=Surv(event_time,event_status)~group,
                              data=sim_data,
                              method="mw",
                              t_star = 0
)
weights_cens<-find_weights(formula=Surv(event_time,event_status)~group,
                              data=sim_data,
                              method="mw",
                              t_star = 0,
                              include_cens=FALSE
)
weights_lr<-find_weights(formula=Surv(event_time,event_status)~group,
                         data=sim_data,
                         method="lr"
)
weights_fh<-find_weights(formula=Surv(event_time,event_status)~group,
                         data=sim_data,
                         method="fh",
                         rho = 0,
                         gamma=1
)

for (weights in 
     list(weights_mw,weights_lr,weights_fh,
          weights_mw_s_star,weights_mw_0,weights_cens)){
  test_that("correct length", {
    expect_true(length(weights)<=2000)
  })
  
  test_that("weights are increasing", {
    expect_true(all(diff(weights)>=0))
  })
  
}

#Test errors

test_that("specifying rho with lr test", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="lr",
                            rho = 0),
               "do not specify rho or gamma for log rank test")
})

test_that("specifying t_star with lr test", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="lr",
                            t_star=4),
               "do not specify t_star or s_star for log rank test ")
})

test_that("specifying rho with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="mw",
                            rho = 0,
                            t_star=4),
               "do not specify rho or gamma for modestly weighted log rank test")
})

test_that("null gamma with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="fh",
                            rho = 0),
               "specify rho and gamma")
})

test_that("specify t_star with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="fh",
                            rho = 0,
                            gamma=1,
                            t_star=4),
               "do not specify t_star or s_star for Fleming-Harrington test")
})
test_that("negative gamma with Fleming Harrington", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="fh",
                            rho = 0,
                            gamma=-1),
               "rho and gamma must be non-negative")
})

test_that("s_star not a probability", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="mw",
                            s_star = 0.5,
                            t_star= 4),
               "must specify either t_star or s_star")
})

test_that("s_star and t_star not a specified", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="mw"),
               "must specify either t_star or s_star")
})

test_that("s_star not a probability", {
  expect_error(find_weights(formula=Surv(event_time,event_status)~group,
                            data=sim_data,
                            method="mw",
                            s_star = -0.5),
               "s_star must a probability")
})

sim_data_strata<-sim_data
sim_data_strata$ecog<-rep(c(0,1),1000)
test_that("don't allow strata",
  {expect_error(find_weights(formula=Surv(event_time,event_status)~group+strata(ecog),
                            data=sim_data_strata,
                            method="mw",
                            s_star = 0.5), "does not account for strata")}
)
