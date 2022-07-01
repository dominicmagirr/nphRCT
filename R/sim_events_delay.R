#' Simulate survival data from a two-arm trial with a delayed separation of survival curves
#'
#' Simulate survival data from a two-arm trial with a delayed separation of survival curves.
#' Survival times on the control arm are simulated from an exponential distribution.
#' Survival times on the experimental arm are simulated from a two-piece exponential distribution.
#'
#' @param{n_c} Number of patients on control treatment.
#' @param{n_e} Number of patients on experimental treatment.
#' @param{delay_e} Length of first period of the experimental arm.
#' @param{lambda_c} Rate of exponential distribution  of the control arm.
#' @param{lambda_e_1} Rate of exponential distribution during first period of the experimental arm.
#' @param{lambda_e_2} Rate of exponential distribution during second period of the experimental arm.
#' @param{rec_period} Parameter used to model recruitment according to power model, see Details.
#' @param{rec_power} Parameter used to model recruitment according to power model, see Details.
#' @param{max_cal_t} Calendar time at which the trial ends, all observations are censored at this time.
#' @return Data frame with columns `event_time`, `event_status` (`1` = event, `0` = censored), and treatment arm indicator `group`.
#' @details
#'
#' Recruitment is modeled using the power model
#' \eqn{P(recruited before T) = (T / rec\_period) ^ rec\_power}, where
#' \eqn{rec_period} is the time at the end of recruitment period, and \eqn{rec_power} controls the rate of recruitment.
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
#' sim_events_delay(
#'   n_c = 5,
#'   n_e = 5,
#'   delay_e = 6,
#'   lambda_c = log(2)/9,
#'   lambda_e_1 = log(2)/9,
#'   lambda_e_2 = log(2)/18,
#'   rec_period = 12,
#'   rec_power = 1,
#'   max_cal_t = 36
#' )
#' @export

sim_events_delay<-function(
  n_c,
  n_e,
  delay_e,
  lambda_c,
  lambda_e_1,
  lambda_e_2,
  rec_period,
  rec_power,
  max_cal_t
){
  rec_c = rec_period * runif(n_c)^(1/rec_power)
  rec_e = rec_period * runif(n_e)^(1/rec_power)
  t_c = rexp(n_c, rate = lambda_c)
  t_1_e = rexp(n_e, rate = lambda_e_1)
  t_2_e = rexp(n_e, rate = lambda_e_2)
  t_e = ifelse(t_1_e < delay_e, t_1_e, delay_e + t_2_e)
  cal_t_c = rec_c + t_c
  cal_t_e = rec_e + t_e

  event_c = +(cal_t_c <= max_cal_t)
  event_e = +(cal_t_e <= max_cal_t)
  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_e)
  df = data.frame(event_time = c(obs_t_c, obs_t_e), event_status = c(event_c,
                                                                     event_e), group = rep(c("control", "experimental"), c(n_c, n_e)))
  df$event_time<-round(df$event_time,2)
  df
}
