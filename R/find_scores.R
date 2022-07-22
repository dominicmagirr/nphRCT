#' Calculate scores
#'
#' Weighted log-rank tests can also be formulated as a permutation test, see Magirr (2021). 
#' This function calculates the scores for the permutation test formulation
#' for different types of weighted log-rank test,
#' the modestly-weighted log-rank test and the Fleming-Harrington (\eqn{\rho},\eqn{\gamma}) 
#' test, in addition to the standard log-rank test.
#'
#' @template formula
#' @template data
#' @param method Character string specifying type of method to calculate scores. Either one of the weighted log-rank tests
#' (log-rank `"lr"`, Fleming-Harrington `"fh"`, modestly weighted `"mw"`) or pseudovalue-based scores (restricted mean survival time `"rmst"`,
#' milestone analysis `"ms"`)
#' @template t_star
#' @template s_star
#' @template rho
#' @template gamma
#' @param tau Parameter \eqn{\tau} in the RMST (`"rmst"`) or milestone analysis  (`"ms"`) test.
#' @return 
#' Data frame. Each row corresponds to an event time (including censoring times if `include_cens=TRUE`).
#' At each time specified in `t_j` the columns indicate
#' - `score_censored` the scores for each censoring at time `t_j` in the permutation test formulation
#' - `score_event` the scores for each event at time `t_j` in the permutation test formulation
#' - `event` the event indicator
#' - `group` the treatment arm indicator
#' - `score` the value of the score for this individual, equal to `score_censored` if the individual is censored 
#' and equal to `score_event` if the individual experienced the event
#' - `standardized_score` the value of `score` standardized to be between -1 and 1
#' - `rank` the rank of the ordered event and censoring times, with higher numbers indicating
#'  later times `t_j`
#' 
#' @details
#'
#' Select which of the tests to perform using argument `method`.
#' For the weighted log-rank tests, the output is calculated as outlined in `vignette("weighted_log_rank_tests", package="wlrt")`.
#' 
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
#' df_scores<-find_scores(formula=Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   method="mw",
#'   t_star = 4
#' )
#' plot(df_scores)
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

find_scores<-function(formula,
                       data,
                       method,
                       t_star = NULL,
                       s_star = NULL,
                       rho = NULL,
                       gamma = NULL,
                       tau=NULL){
  method<-match.arg(method,c("mw","lr","fh","rmst","ms"))
  
  #weighted log rank tests
  if(method %in% c("mw","lr","fh")){
    w<-find_weights(formula=formula,
                   data=data,
                   method=method,
                   t_star = t_star,
                   s_star = s_star,
                   rho = rho,
                   gamma = gamma,
                   include_cens=TRUE)
    
    df <- find_at_risk(formula=formula,
                               data=data,
                               include_cens=TRUE)
    df$score_censored<-with(df,-cumsum(w*(n_event/n_risk)))
    df$score_event<-with(df,score_censored+w)
    
    formula_vars <- all.vars(formula)
    time_col <- formula_vars[1]
    status_col <- formula_vars[2]
    group_col<-formula_vars[-(1:2)]
    
    data[["t_j"]]<-data[[time_col]]
    data[["event"]]<-data[[status_col]]
    data[["group"]]<-data[[group_col]]
    
    df<-merge(x = df[,c("t_j","score_event","score_censored")], y = data[,c("t_j","event","group")], by = "t_j", all.x = TRUE)

    df$score <- with(df, event * score_event + (1-event) * score_censored)
    df$score_censored<-NULL
    df$score_event<-NULL
  }
  
  #Restricted mean survival time
  if(method=="rmst"){
    df <-find_pseudovalues(formula=formula,
                      data=data,
                      method = "rmst",
                      tau = tau)
  }
  
  #Milestone
  if(method=="ms"){
    df <-find_pseudovalues(formula=formula,
                                data=data,
                                method = "ms",
                                tau = tau)
  }

  ###########################################
  ## pick out score (event or censored) (DM)
  ## standardize scores to (-1, 1) (DM)
  max_a <- max(df$score)
  min_a <- min(df$score)
  A = 2 / (max_a - min_a)
  B = 1 - A * max_a
  df$standardized_score <- df$score * A + B
  ############################################
  
  # df$rank<-factor(df$t_j)
  # levels(df$rank)<-paste0("(",1:length(unique(df$t_j)),")")
  out<-list(df=df)
  class(out)<-"df_score"
  out
}
args<-c()
#' @export
plot.df_score<-function(x,...){
  df<-x[["df"]]
  df[["group"]]<-as.factor(df[["group"]])
  gl1=levels(df[["group"]])[1]
  gl2=levels(df[["group"]])[2]
  df$event<-as.factor(df[["event"]])

  args <- list(col="Arm",x="Time",y="Score")
  inargs <- list(...)
  args[names(inargs)] <- inargs
  labels<-do.call(ggplot2::labs,args)
  
  means<-data.frame(intercept=c(mean(df[df[["group"]]==gl1,"standardized_score"]),mean(df[df[["group"]]==gl2,"standardized_score"])),
                         group=c(gl1,gl2))
  ggplot2::ggplot(df, ggplot2::aes_string(x="t_j", y="standardized_score",col="group")) + ggplot2::geom_point(alpha=0.3) +
    ggplot2::ylim(-1.1,1.1)+labels+ 
    ggplot2::geom_hline(ggplot2::aes_string(yintercept="intercept",col="group"),linetype="dashed",data=means)
  }
