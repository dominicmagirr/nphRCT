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
  n_c=5000,
  n_e=5000,
  max_cal_t = 36
)

#with censoring included
out_cens<-find_at_risk(formula=Surv(event_time,event_status)~group,
             data=sim_data,
             include_cens=TRUE)
#with censoring excluded
out_no_cens<-find_at_risk(formula=Surv(event_time,event_status)~group,
             data=sim_data,
             include_cens=FALSE)

for (out in list(out_cens,out_no_cens)){
  test_that("check time is strictly increasing", {
    expect_true(all(diff(out$t_j)>0))
  })
  
  test_that("check number at risk is decreasing", {
    out_n_risk<-out[,grep("n_risk",names(out))]
    test_decreasing<-function(x){all(diff(x)<=0)}
    expect_true(all(apply(out_n_risk,2,test_decreasing)))
  })  
  
  test_that("check events is greater than / equal to zero", {
    out_n_event<-out[,grep("n_event",names(out))]
    test_greater_zero<-function(x){is.numeric(x) && all(x>=0)}
    expect_true(all(apply(out_n_event,2,test_greater_zero)))
  })
}



#Test errors
sim_data_strata<-sim_data
sim_data_strata$ecog<-rep(c(0,1),1000)

test_that("null gamma with Fleming Harrington", {
  expect_error(find_at_risk(formula=Surv(event_time,event_status)~group+strata(ecog),
                            data=sim_data_strata),
               "Function does not account for strata")
})

sim_data_groups<-sim_data
sim_data_groups$group<-rep(1:4,each=2500)
test_that("negative gamma with Fleming Harrington", {
  expect_error(find_at_risk(formula=Surv(event_time,event_status)~group,
                            data=sim_data_groups),
               "Only 2 treatment groups")
})
