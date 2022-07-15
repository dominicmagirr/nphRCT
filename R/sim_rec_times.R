#' Simulate recruitment times
#'
#' Simulate recruitment times from either a power model or a piecewise constant model
#'
#' @param{rec_model} Character string specifying the type of recruitment model. Either the power model `"power"` 
#' or piecewise constant `"pw_constant"`.
#' @param{n} Number of recruitment times to simulate
#' @param{rec_power} Parameter used to model recruitment according to power model, see Details. 
#' @param{rec_period} Parameter used to model recruitment according to power model, see Details. 
#' @param{rec_rate}  Parameter used to model recruitment according to piecewise constant model, see Details. 
#' @param{rec_duration} Parameter used to model recruitment according to piecewise constant model, see Details. 
#' @return 
#' @details
#'
#' Recruitment is modeled using either the power model or the piecewise constant model. 
#' 
#' The power model is defined as:
#' \eqn{P(recruited before T) = (T / rec_period) ^ rec_power}, where
#' \eqn{rec_period} is the time at the end of recruitment period, and \eqn{rec_power} controls the rate of recruitment.
#'
#' Alternatively, recruitment can be modelled using the piecewise constant model. 
#' In the simple case with only one time period defined in `rec_duration`, the times between each of the individuals entering follow-up
#' are samples from the exponential distribution with rate parameter \eqn{\lambda},
#' \eqn{f(t)=\lambda exp(-\lambda t)}. The function returns `n` recruitment times, regardless of the 
#' length of duration `rec_duration.`
#' 
#' In the case with multiple time periods defined in `rec_duration`, the number of events in each period is sampled from the 
#' Poisson distribution \eqn{P(K=k)=\lambda^k \exp{-\lambda}/k!}, where \eqn{k} is the number of events. The rate parameter 
#' \eqn{\lambda} is equal to `rec_rate` multiplied by the duration of the time period in `rec_duration`. The \eqn{K} recruitment times 
#' are then sampled uniformly from the corresponding time period.
#' 
#' In the case with multiple time periods, fixing the number of events `n` is optional. In the case that `n` is fixed, and 
#' insufficient recruitment times have been simulated by the end of the last time period, the additional recruitment times will be 
#' simulated after the end of the last time period using the method described above for when there is one time period.
#'  
#'  
#' @examples
#' library(wlrt)
#' set.seed(1)
#' #power model
#' sim_rec_times(n=10,
#'   rec_model= "power",
#'   rec_period=20,
#'   rec_power=1)
#'
#' #piecewise constant model, where duration is set
#' sim_rec_times(rec_model="pw_constant",
#'           rec_rate=c(1,2),
#'           rec_duration=c(12,8))
#' #piecewise constant model, where number of events is set
#' sim_rec_times(rec_model="pw_constant",
#'           rec_rate=c(1,2),
#'           rec_duration=c(12,8),
#'           n=10)
#'    
#' @export


sim_rec_times<-function(rec_model,
                        n=NULL,
                        rec_power=NULL,
                        rec_period=NULL,
                        rec_rate=NULL,
                        rec_duration=NULL){
  rec_model <- match.arg(rec_model, c("power", "pw_constant"))
  
  if(rec_model=="power"){
    rec <- rec_period * runif(n)^(1/rec_power)
    return(rec)
  }
  if(rec_model=="pw_constant"){
    if(any(rec_rate<0)){stop("rec_rate should be non-negative")}
    if(length(rec_rate)==1){#simple case with only one rate
      rec<-cumsum(stats::rexp(n=n,rate=rec_rate))
      return(rec)
    }else{#piecewise
      if(length(rec_duration)!=length(rec_rate)){stop("Lengths of rec_duration and rec_lambda should match")}
      n_periods<-length(rec_duration)
      df<-data.frame(rate=rec_rate,
                     duration=rec_duration,
                     period=1:n_periods,
                     finish=cumsum(rec_duration),
                     lambda=rec_duration*rec_rate,
                     origin=c(0,rec_duration[-n_periods]))
      df$N<-sapply(df$lambda,function(x){stats::rpois(n=1,lambda = x)})
      if (sum(df$N)==0){
        if (is.null(n)){return(NULL)}
        if (df$rate[n_periods]==0) stop("Please specify positive rec_lambda for the last period; otherwise enrollment cannot finish.")
        rec<-c(cumsum(stats::rexp(n,rate=df$rate[n_periods]))+df$finish[n_periods])
        return(rec)
      }else{
        rec<-unlist(apply(df,1,function(x){sort(stats::runif(n=x[["N"]],min=x[["origin"]],max=x[["finish"]]))}))
        if (is.null(n)){return(rec)} # if n not specified, return generated times
        if (length(rec) >= n){return(rec[1:n])} # if n already achieved, return first n observations
        # stop with error message if enrollment has not finished but enrollment rate for last period is less or equal with 0
        if (df$rate[n_periods]==0){stop("Please specify positive rec_lambda for the last period; otherwise enrollment cannot finish.")}
        # Otherwise, return inter-arrival exponential times
        rec<-c(rec, cumsum(stats::rexp(n-nrow(rec),rate=df$rate[n_periods]))+df$finish[n_periods])
        return(rec)
      }
    }
  }
}