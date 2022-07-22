#' Simulate survival data from a two-arm trial with a delayed separation of survival curves
#'
#' Simulate survival data from a two-arm trial with a delayed separation of survival curves.
#' Survival times on the control arm are simulated from an exponential distribution.
#' Survival times on the experimental arm are simulated from a two-piece exponential distribution.
#'
#' @param delay_e Length of first period of the experimental arm.
#' @param lambda_c Rate of exponential distribution  of the control arm.
#' @param lambda_e_1 Rate of exponential distribution during first period of the experimental arm.
#' @param lambda_e_2 Rate of exponential distribution during second period of the experimental arm.
#' @param rec_times_c Recruitment times on the control arm with each value corresponding to one individual
#' @param rec_times_e Recruitment times on the experimental arm with each value corresponding to one individual 
#' @param max_cal_t Calendar time at which the trial ends, all observations are censored at this time.
#' @return Data frame with columns `event_time`, `event_status` (`1` = event, `0` = censored), and treatment arm indicator `group`.
#' @details
#'
#' Survival times are simulated from an exponential distribution with rate parameter \eqn{\lambda},
#' \eqn{f(t)=\lambda exp(-\lambda t)}. This distribution has a median value of \eqn{log(2)/\lambda};
#' this can be a useful fact when setting the rates `lambda_c`, `lambda_e_1`, `lambda_e_2`.
#'
#' All observations are censored at the calendar time defined in argument `max_cal_t`.
#'
#' @examples
#' library(wlrt)
#' set.seed(1)
#' rec_c <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=5)
#' rec_e <- sim_rec_times(rec_model="power",rec_period=12,rec_power=1,n=5)
#' sim_data <- sim_events_delay(
#'   delay_e = 6,
#'   lambda_c = log(2)/9,
#'   lambda_e_1 = log(2)/9,
#'   lambda_e_2 = log(2)/18,
#'   rec_times_c = rec_c,
#'   rec_times_e = rec_e,
#'   max_cal_t = 36
#' )
#' @export



sim_events_delay<-function(
  delay_e,
  lambda_c,
  lambda_e_1,
  lambda_e_2,
  rec_times_c,
  rec_times_e,
  max_cal_t
){
  if(any(c(rec_times_c,rec_times_e) > max_cal_t)|any(c(rec_times_c,rec_times_e) < 0)){stop("patients must enroll between 0 and max_cal_t")}

  n_c=length(rec_times_c)
  n_e=length(rec_times_e)
  t_c = stats::rexp(n_c, rate = lambda_c)
  
  t_1_e = stats::rexp(n_e, rate = lambda_e_1)
  t_2_e = stats::rexp(n_e, rate = lambda_e_2)
  t_e = ifelse(t_1_e < delay_e, t_1_e, delay_e + t_2_e)
  cal_t_c = rec_times_c + t_c
  cal_t_e = rec_times_e + t_e

  event_c = +(cal_t_c <= max_cal_t)
  event_e = +(cal_t_e <= max_cal_t)
  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_times_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_times_e)
  df = data.frame(event_time = c(obs_t_c, obs_t_e), event_status = c(event_c,
                                                                     event_e), group = rep(c("control", "experimental"), c(n_c, n_e)))
  df$event_time<-round(df$event_time,2)
  df
}

