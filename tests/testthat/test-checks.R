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
sim_data<-sim_data[sim_data$event_time!=0,]
sim_data$event_start<-0
test_that("formula is correct", {
  expect_error(check_formula(formula=Surv(time=event_start,time2=event_time,event_status)~group,
                            data=sim_data),
               "Censoring type should be right censoring")
  expect_error(check_formula(formula=survival::Surv2(event_time,event_status)~group,data=sim_data),
               "formula has an incorrect format")
})


sim_data$event_start<-NULL
sim_data_0 <- sim_data
sim_data_0$ecog=0
sim_data_1 <- sim_events_delay(
  event_model=list(
    duration_c = 36,
    duration_e = c(6,30),
    lambda_c = log(2)/6,
    lambda_e = c(log(2)/6,log(2)/12)
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
sim_data_1$ecog=1
sim_data_strata<-rbind(sim_data_0,sim_data_1)

sim_data_NA_1<-sim_data_NA_2<-sim_data_NA_3<-sim_data_NA_4<-sim_data_strata
sim_data_NA_1[1,"event_time"]<-NA
sim_data_NA_2[1,"event_status"]<-NA
sim_data_NA_3[1,"group"]<-NA
sim_data_NA_4[1,"ecog"]<-NA

for (sim_data_NA in list(sim_data_NA_1,sim_data_NA_2,sim_data_NA_3,sim_data_NA_4)){
  test_that("does not work with missing data", {
    expect_error(check_formula(formula=Surv(event_time,event_status)~group,
                               data=sim_data_NA_1),
                 "NAs in data set")
  })
}

