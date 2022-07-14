#' Calculate weights
#'
#' This function can perform two types of weighted log-rank test,
#' the modestly-weighted log-rank test and the Fleming-Harrington (\eqn{\rho},\eqn{\gamma}) test, in addition to the standard log-rank test.
#'
#' @template formula
#' @template data
#' @template wlr
#' @template t_star
#' @template s_star
#' @template rho
#' @template gamma
#' @template include_cens
#' @return Vector of weights in the weighted log-rank test.
#' The weights correspond to the ordered, distinct event times (and censoring times if 
#' `include_cens=TRUE`).
#'
#' @details
#'
#' Select which of the three tests to perform using argument `wlr`.
#' The output is calculated as outlined in `vignette("weighted_log_rank_tests", package="wlrt")`.
#' 
#' For the weighted log rank tests, the weights are
#' multiplied by the observed minus expected value at each event time.
#' The weights that correspond to censoring times are therefore not used in the weighted log-rank test, 
#' however they are used in the `find_scores` function.
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
#' #example setting t_star
#' find_weights(formula=Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   wlr="mw",
#'   t_star = 4
#' )
#' @references
#' Magirr, D. (2021).
#' Non-proportional hazards in immuno-oncology: Is an old perspective needed?.
#' Pharmaceutical Statistics, 20(3), 512-527.
#' DOI: 10.1002/pst.2091
#'
#' Magirr, D. and Burman, C.F., 2019.
#' Modestly weighted logrank tests.
#' Statistics in medicine, 38(20), 3782-3790.
#' @export

find_weights<-function(formula,
                           data,
                           wlr,
                           t_star = NULL,
                           s_star = NULL,
                           rho = NULL,
                           gamma = NULL,
                           include_cens=FALSE){

  check_formula(formula=formula,data=data)
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  Terms <- terms(formula,"strata")
  strata_index <- survival::untangle.specials(Terms,"strata")$terms
  if(length(strata_index)>0){stop("Function does not account for strata")}
  
  #Find weights
  wlr <- match.arg(wlr, c("lr", "fh", "mw"))
  Surv <- survival::Surv
  formula_km<-as.formula(paste0("Surv(",time_col,",",status_col,")~1"))
  km_fit <- survival::survfit(formula_km,
                                data = data,
                                timefix = FALSE)
  if(include_cens==TRUE){
    t_j <- km_fit$time
    S_hat <- km_fit$surv
    n_risk <- km_fit$n.risk
    n_event <- km_fit$n.event
    n_censor <- km_fit$n.censor
  }else{
    events<-km_fit$n.event > 0
    t_j <- km_fit$time[events]
    S_hat <- km_fit$surv[events]
    n_risk <- km_fit$n.risk[events]
    n_event <- km_fit$n.event[events]
    n_censor <- km_fit$n.censor[events]
  }
  S_hat_minus <- c(1, S_hat[1:(length(S_hat) - 1)])
  
  if (wlr == "lr") {
    check_lr(rho=rho,gamma=gamma,t_star=t_star,s_star=s_star)
    w <- rep(1, length(S_hat_minus))
  }
  if (wlr == "fh") {
    check_fh(rho=rho,gamma=gamma,t_star=t_star,s_star=s_star)
    w <- S_hat_minus ^ rho * (1 - S_hat_minus) ^ gamma
  }
  if (wlr == "mw") {
    check_mw(rho=rho,gamma=gamma,t_star=t_star,s_star=s_star)
    if (!is.null(t_star)){
      if (any(t_j < t_star)) {
        w <- pmin(1 / S_hat_minus,
                  1 / S_hat[max(which(t_j < t_star))])
      }
      else {
        w <- rep(1, length(S_hat_minus))
      }
    }
    else {
      if(s_star<0 | s_star>1){stop("s_star must a probability (between 0 and 1)")}
      w <- pmin(1 / S_hat_minus,
                1 / s_star)
    }
  }

  w
}
