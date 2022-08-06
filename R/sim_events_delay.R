#' Simulate survival data from a two-arm trial
#'
#' Simulate survival data from a two-arm trial with
#' survival times on the control arm and experimental arm simulated from an exponential distribution 
#' or piecewise exponential distribution.
#' @param event_model List containing information to simulate event times, including:
#' - `duration_c` Vector of durations corresponding to each of the periods of the control arm.
#' - `duration_e` Vector of durations corresponding to each of the periods of the experimental arm.
#' - `lambda_c` Vector of parameters \eqn{\lambda} in the exponential distribution corresponding to each of the periods of the control arm.
#' - `lambda_e` Vector of parameters \eqn{\lambda} in the exponential distribution corresponding to each of the periods of the experimental arm.
#' @param recruitment_model List containing information to simulate recruitment times, including:
#' - `rec_model` Character string specifying the type of recruitment model. Either the power model `"power"` 
#' or piecewise constant model `"pw_constant"`.
#' - `rec_power` Parameter used to model recruitment according to power model, see Details. 
#' - `rec_period` Parameter used to model recruitment according to power model, see Details. 
#' - `rec_rate` Parameter used to model recruitment according to piecewise constant model, see Details. 
#' - `rec_duration` Parameter used to model recruitment according to piecewise constant model, see Details. 
#' @param max_cal_t Calendar time at which the trial ends, all observations are censored at this time.
#' @param n_c Number of individuals on the control arm 
#' @param n_e Number of individuals on the event arm 
#' @return Data frame with columns `event_time`, `event_status` (`1` = event, `0` = censored), and treatment arm indicator `group`.
#' @details
#'
#' Survival times are simulated from an exponential distribution with rate parameter \eqn{\lambda},
#' \eqn{f(t)=\lambda exp(-\lambda t)}. This distribution has a median value of \eqn{log(2)/\lambda};
#' this can be a useful fact when setting the rates `lambda_c` and `lambda_e`. 
#' The survival times can be simulated from a piecewise exponential distribution, setting one/multiple
#' durations and \eqn{\lambda} parameters for the control and experimental arms.
#' 
#' Recruitment is modeled using either the power model or the piecewise constant model. 
#' 
#' The power model is defined as:
#' \eqn{P(recruited\_before\_T) = (T / rec\_period) ^ {rec\_power}}, where
#' \eqn{rec\_period} is the time at the end of recruitment period, and \eqn{rec\_power} controls the rate of recruitment.
#' 
#' Alternatively, recruitment can be modelled using the piecewise constant model. 
#' In the simple case with only one time period defined in `rec_duration`, the times between each of the individuals entering follow-up
#' are samples from the exponential distribution with rate parameter \eqn{\lambda},
#' \eqn{f(t)=\lambda exp(-\lambda t)}. The number of recruitment times defined in `n_c` or `n_e` is returned, regardless of the 
#' length of duration `rec_duration.`
#' 
#' In the case with multiple time periods defined in `rec_duration`, the number of events in each period is sampled from the 
#' Poisson distribution \eqn{P(K=k)=\lambda^k \exp{(-\lambda}/k!)}, where \eqn{k} is the number of events. The rate parameter 
#' \eqn{\lambda} is equal to `rec_rate` multiplied by the duration of the time period in `rec_duration`. The recruitment times 
#' are then sampled uniformly from the corresponding time period. In the case that
#' insufficient recruitment times have been simulated by the end of the last time period, the additional recruitment times will be 
#' simulated after the end of the last time period.
#'  
#' All observations are censored at the calendar time defined in argument `max_cal_t`.
#'
#' @examples
#' library(nphRCT)
#' set.seed(1)
#' sim_data <- sim_events_delay(
#'   event_model=list(
#'     duration_c = 36,
#'     duration_e = c(6,30),
#'     lambda_c = log(2)/9,
#'     lambda_e = c(log(2)/9,log(2)/18)
#'   ),
#'   recruitment_model=list(
#'     rec_model="power",
#'     rec_period = 12,
#'     rec_power = 1
#'   ),
#'   n_c=50,
#'   n_e=50,
#'   max_cal_t = 36
#' )
#' @export

sim_events_delay<-function(
  event_model,
  recruitment_model,
  n_c,
  n_e,
  max_cal_t
){
  rec_model <- match.arg(recruitment_model$rec_model, c("power", "pw_constant"))
  n_total<-n_c+n_e
  if(rec_model=="power"){
    rec_period<-recruitment_model$rec_period
    rec_power<-recruitment_model$rec_power
    rec <- rec_period * stats::runif(n_total)^(1/rec_power)
  }
  if(rec_model=="pw_constant"){
    rec_rate<-recruitment_model$rec_rate
    if(any(rec_rate<0)){stop("rec_rate should be non-negative")}
    if(length(rec_rate)==1){#simple case with only one rate
      rec<-cumsum(stats::rexp(n=n_total,rate=rec_rate))
    }else{#piecewise
      rec_duration<-recruitment_model$rec_duration
      if(length(rec_duration)!=length(rec_rate)){stop("Lengths of rec_duration and rec_rate should match")}
      n_periods<-length(rec_duration)
      df<-data.frame(rate=rec_rate,
                     duration=rec_duration,
                     period=1:n_periods,
                     finish=cumsum(rec_duration),
                     lambda=rec_duration*rec_rate,
                     origin=c(0,cumsum(rec_duration)[-n_periods]))
      df$N<-sapply(df$lambda,function(x){stats::rpois(n=1,lambda = x)})
      if (sum(df$N)==0){
        if (df$rate[n_periods]==0) stop("Please specify positive rec_rate for the last period; otherwise enrollment cannot finish.")
          rec<-c(cumsum(stats::rexp(n_total,rate=df$rate[n_periods]))+df$finish[n_periods])
        }else{
          rec<-unlist(apply(df,1,function(x){sort(stats::runif(n=x[["N"]],min=x[["origin"]],max=x[["finish"]]))}))
          if (length(rec) >= n_total){rec<-rec[1:n_total]} # if n already achieved, return first n observations
          # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
          else{if (df$rate[n_periods]==0){stop("Please specify positive rec_rate for the last period; otherwise enrollment cannot finish.")}
          # Otherwise, return inter-arrival exponential times
          rec<-c(rec, cumsum(stats::rexp(n_total-nrow(rec),rate=df$rate[n_periods]))+df$finish[n_periods])
          }
        }
      }
    }
  
  sample_n_c <- sample(1:n_total, n_c, replace=FALSE)
  sample_n_e <- sample(setdiff(1:n_total, sample_n_c), n_e, replace=FALSE)
  rec_times_c<-rec[sample_n_c]
  rec_times_e<-rec[sample_n_e]
  
  if(any(c(rec_times_c,rec_times_e) > max_cal_t)|any(c(rec_times_c,rec_times_e) < 0)){stop("patients must enroll between 0 and max_cal_t")}
  
  duration_c<-event_model$duration_c
  duration_e<-event_model$duration_e
  lambda_c<-event_model$lambda_c
  lambda_e<-event_model$lambda_e
  if(sum(duration_c)<max_cal_t | sum(duration_e)<max_cal_t){stop("The total length of durations in duration_c and duration_e must not be less than max_cal_t")}
  #control
  t_c <- t_c_period <- sapply(lambda_c,function(lambda){stats::rexp(n_c, rate = lambda)})
  duration_c[length(duration_c)] <- Inf
  t_c_period<-apply(t_c_period,1, function(i){min(which(i<duration_c))})
  t_c<-sapply(seq_along(duration_c), function(i){ifelse(t_c[,i]>duration_c[i],duration_c[i],t_c[,i])})
  t_c<-sapply(1:n_c, function(i){sum(t_c[i,1:t_c_period[i]])})
  #experimental
  t_e = t_e_period = sapply(lambda_e,function(lambda){stats::rexp(n_e, rate = lambda)})
  duration_e[length(duration_e)] = Inf
  t_e_period=apply(t_e_period,1, function(i){min(which(i<duration_e))})
  t_e<-sapply(seq_along(duration_e), function(i){ifelse(t_e[,i]>duration_e[i],duration_e[i],t_e[,i])})
  t_e<-sapply(1:n_e, function(i){sum(t_e[i,1:t_e_period[i]])})
  
  cal_t_c = rec_times_c + t_c
  cal_t_e = rec_times_e + t_e
  
  event_c = +(cal_t_c <= max_cal_t)
  event_e = +(cal_t_e <= max_cal_t)
  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_times_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_times_e)
  df = data.frame(event_time = c(obs_t_c, obs_t_e), 
                  event_status = c(event_c,event_e), 
                  group = c(rep("control",n_c), rep("experimental",n_e)))
  df$event_time<-round(df$event_time,2)
  df
}

