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
#' @template wlr
#' @template t_star
#' @template s_star
#' @template rho
#' @template gamma
#' @return Data frame containing output from function `find_weights` and additional columns:
#' 
#' - `score_cens` the scores for each censoring at time `t_j` in the permutation test formulation
#' - `score_event` the scores for each event at time `t_j` in the permutation test formulation
#' - `rank` the rank of the ordered event and censoring times, the higher the number the earlier the time `t_j`
#' @details
#'
#' Select which of the three tests to perform using argument `wlr`.
#' The output is calculated as outlined in `vignette("weighted_log_rank_tests", package="wlrt")`.
#'
#' @examples
#' library(wlrt)
#' set.seed(1)
#' sim_data <- sim_events_delay(
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
#' find_scores(formula=Surv(event_time,event_status)~group,
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

find_scores<-function(formula,
                       data,
                       wlr,
                       t_star = NULL,
                       s_star = NULL,
                       rho = NULL,
                       gamma = NULL,
                       ...){
  df<-find_weights(data=data,
                 formula=formula,
                 wlr=wlr,
                 t_star = t_star,
                 s_star = s_star,
                 rho = rho,
                 gamma = gamma,
                 include_cens=TRUE)

  df$score_cens<-with(df,-cumsum(w*(n_event/n_risk)))
  df$score_event<-with(df,score_cens+w)
  df$rank<-paste0("(",nrow(df):1,")")

  out<-list(df=df)
  class(out)<-"wlrt_score"
  out
}

#' @export
plot.wlrt_score<-function(x,...){
  df<-x$df
  df$x_pos<-1:nrow(df)
  df_cens<-df[df$n_censor>0,]
  df_event<-df[df$n_event>0,]
  
  args <- list(ylim=c(min(df_cens$score_cens, df_event$score_event)-0.5,
                      max(df_cens$score_cens, df_event$score_event)+0.5))
  inargs <- list(...)
  args[names(inargs)] <- inargs
  
  do.call(plot,c(list(x=df$x_pos,
                      type="n",xlab="Rank",xaxt="n",ylab="Score"),args))
  
  axis(side=1, labels=df$rank, at=df$x_pos)
  points(x=df_cens$x_pos,y=df_cens$score_cens,pch=21,bg="white")
  points(x=df_event$x_pos,y=df_event$score_event,pch=21,bg="black")
  
}
