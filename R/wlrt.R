#' Weighted log-rank test
#'
#' This function can perform two types of weighted log-rank test,
#' the modestly-weighted log-rank test and the Fleming-Harrington (\eqn{\rho},\eqn{\gamma}) test, in addition to the standard log-rank test.
#'
#' @template formula
#' @template data
#' @param method Character string specifying type of weighted log-rank test. 
#' Either `"lr"` for a standard log-rank test, `"mw"` for a modestly-weighted log-rank test, 
#' or `"fh"` for the Fleming-Harrington rho-gamma family.
#' @template t_star
#' @template s_star
#' @template rho
#' @template gamma
#' @param timefix Deal with floating point issues (as in the survival package). Default is TRUE. May need to set FALSE for simulated data.
#' @return List containing the outcome of the weighted log-rank test.
#' - `u` is the test statistic U for the weighted log-rank test
#' - `v_u` is the variance of test statistic U
#' - `z` is the Z-score
#' - `trt_group` indicates which of the treatment arms the test statistic U corresponds to
#'
#' In the presence of multiple strata, the results of the test on each individual strata is returned, in addition to the combined test that was proposed by 
#' Magirr and Jiménez (2022), see `vignette("weighted_log_rank_tests", package="wlrt")`.
#' @details
#'
#' Select which of the three tests to perform using argument `method`.
#' The output is calculated as outlined in `vignette("weighted_log_rank_tests", package="wlrt")`.
#'
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
#' #example setting t_star
#' wlrt(formula=Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   method="mw",
#'   t_star = 4
#' )
#' #example setting s_star
#' wlrt(formula=Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   method="mw",
#'   s_star = 0.5
#' )
#' #example with 1 strata
#' sim_data_0 <- sim_data
#' sim_data_0$ecog=0
#' sim_data_1 <- sim_events_delay(
#'   event_model=list(
#'     duration_c = 36,
#'     duration_e = c(6,30),
#'     lambda_c = log(2)/6,
#'     lambda_e = c(log(2)/6,log(2)/12)
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
#' sim_data_1$ecog=1
#' sim_data_strata<-rbind(sim_data_0,sim_data_1)
#' wlrt(formula=Surv(event_time,event_status)~group+strata(ecog),
#'   data=sim_data_strata,
#'   method="mw",
#'   t_star = 4
#' )
#' #example with 2 strata
#' sim_data_strata_2<-cbind(sim_data_strata,sex=rep(c("M","F"),times=100))
#' wlrt(formula=Surv(event_time,event_status)~group+strata(ecog)+strata(sex),
#'   data=sim_data_strata_2,
#'   method="mw",
#'   t_star = 4
#' )
#' @references
#' Magirr, D. (2021).
#' Non-proportional hazards in immuno-oncology: Is an old perspective needed?.
#' Pharmaceutical Statistics, 20(3), 512-527.
#' <doi:10.1002/pst.2091>
#'
#' Magirr, D. and Burman, C.F., 2019.
#' Modestly weighted logrank tests.
#' Statistics in medicine, 38(20), 3782-3790.
#' 
#' Magirr, D. and Jiménez, J. (2022)
#' Stratified modestly-weighted log-rank tests in settings with an anticipated delayed separation of survival curves
#' PREPRINT at <https://arxiv.org/abs/2201.10445>
#' @export

wlrt <- function(formula,
                 data,
                 method,
                 t_star = NULL,
                 s_star = NULL,
                 rho = NULL,
                 gamma = NULL,
                 timefix = TRUE) {
  
  check_formula(formula=formula,data=data)
  Surv<-survival::Surv
  
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  terms_vars<-formula_vars[-(1:2)]
  Terms <- stats::terms(formula,"strata")
  strata_index <- survival::untangle.specials(Terms,"strata")$terms

  if (length(strata_index)>0){
    strata_col<-terms_vars[strata_index]
    group_col<-terms_vars[-strata_index]  
  }else{
    group_col<-terms_vars
  }
  
  if (length(group_col)!=1) {
     stop("Formula must contain only one treatment arm indicator. All other terms must be specified as strata.")
  }
  data[[group_col]]<-as.factor(data[[group_col]])
    
  if (length(strata_index)==0) {
    
    out<-wlrt_strata(formula=formula,
           data=data,
           method=method,
           t_star = t_star,
           s_star = s_star,
           rho = rho,
           gamma = gamma,
           timefix = timefix)
    
  }
  else{
    formula<-stats::as.formula(paste0("Surv(",time_col,",",status_col,") ~ ",group_col))
    for (str in strata_col){data[[str]]<-paste0(str, data[[str]])}
    data_strata<-split(data, data[,strata_col])
    if (min(purrr::map_dbl(data_strata, function(x)dim(x)[1])) < 5){
      stop("Minimum stratum size is 5. Consider merging strata.")
    }
    test_w<-purrr::map_df(data_strata,
                          wlrt_strata,
                          formula=formula,
                          method=method,
                          t_star = t_star,
                          s_star = s_star,
                          rho = rho,
                          gamma = gamma,
                          timefix = timefix)
    
    test_lr<-purrr::map_df(data_strata,
                           wlrt_strata,
                           formula=formula,
                           method="lr",
                           timefix = timefix)
    
    u_tilde_w_z<-sum(sqrt(test_lr$v_u)*test_w$z)
    v_tilde_w_z<-sum(test_lr$v_u)
    z_tilde_w_z<-u_tilde_w_z/sqrt(v_tilde_w_z)
    
    out <- list(by_strata = data.frame(strata=names(data_strata),test_w),
                combined = data.frame(
                  u = u_tilde_w_z,
                  v = v_tilde_w_z,
                  z = z_tilde_w_z,
                  trt_group = sort(unique(data[[group_col]]))[2]
                ))
  }
  out
  
}


wlrt_strata <- function(formula,
                        data,
                        method,
                        t_star = NULL,
                        s_star = NULL,
                        rho = NULL,
                        gamma = NULL,
                        timefix = TRUE) {
  
  formula_vars <- all.vars(formula)
  terms_vars<-formula_vars[-(1:2)]
  
  #calculate weights
  w<-find_weights(formula=formula,
                  data=data,
                  method=method,
                  t_star = t_star,
                  s_star = s_star,
                  rho = rho,
                  gamma = gamma,
                  include_cens=FALSE,
                  timefix = timefix)

  #at risk tables
  df_events<-find_at_risk(formula=formula,
                          data=data,
                          include_cens=FALSE,
                          timefix = timefix)
  
  n_event_g <- as.matrix(df_events[, grep("n_event_",names(df_events))])
  n_event<-df_events$n_event
  n_risk_g <- as.matrix(df_events[, grep("n_risk_",names(df_events))])
  n_risk<-df_events$n_risk
  
  # test statistics
  observed <- w %*% n_event_g
  expected <- w %*% (n_risk_g * (n_event / n_risk))
  u <- (observed - expected)[, 2]
  
  v_u <-w ^ 2 * (apply(n_risk_g, 1, prod) * n_event * (n_risk - n_event)) / (n_risk ^ 2 * (n_risk - 1))
  v_u[n_risk == 1] <- 0
  v_u<-sum(v_u)
  
  z = u / sqrt(v_u)
  
  out <- data.frame(
    u = as.vector((observed - expected)[, 2]),
    v_u = as.vector(v_u),
    z = as.vector(z),
    trt_group = sort(unique(data[[terms_vars]]))[2]
  )
  out
}
