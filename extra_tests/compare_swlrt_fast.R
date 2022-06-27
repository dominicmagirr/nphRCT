source("/home/barreis3/Documents/Gitlab/stratified_weighted_paper/wlrt_fast.R")
source("/home/barreis3/Documents/Gitlab/stratified_weighted_paper/swlrt_fast.R")
devtools::load_all()
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
#example with strata
sim_data_0 <- sim_data
sim_data_0$ecog=0
sim_data_1 <- sim_events_delay(
  n_c = 5,
  n_e = 5,
  delay_e = 6,
  lambda_c = log(2)/6,
  lambda_e_1 = log(2)/6,
  lambda_e_2 = log(2)/12,
  rec_period = 12,
  rec_power = 1,
  max_cal_t = 36
)
sim_data_1$ecog=1
sim_data_strata<-rbind(sim_data_0,sim_data_1)
profvis::profvis(
swlrt_fast(df=sim_data_strata,
            trt_colname="group",
            time_colname="event_time",
            event_colname="event_status",
            strat_colname = "ecog",
            wlr = "mw",
            t_star = 4))
profvis::profvis(
wlrt(formula=Surv(event_time,event_status)~group+strata(ecog),
     data=sim_data_strata,
     wlr="mw",
     t_star = 4
))
