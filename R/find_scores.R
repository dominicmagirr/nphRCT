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
#' @return 
#' Data frame. Each row corresponds to an event time (including censoring times if `include_cens=TRUE`).
#' At each time specified in `t_j` the columns indicate
#' - `score_censored` the scores for each censoring at time `t_j` in the permutation test formulation
#' - `score_event` the scores for each event at time `t_j` in the permutation test formulation
#' - `event` the event indicator
#' - `group` the treatment arm indicator
#' - `score` the value of the score for this individual, equal to `score_censored` if the individual is censored 
#' and equal to `score_event` if the individual experienced the event
#' - `standardized_score` the value of `score` standardised t obe between -1 and 1
#' - `rank` the rank of the ordered event and censoring times, with higher numbers indicating
#'  later times `t_j`
#'
#' Data frame containing output from function `find_weights` and also additional columns:
#' 
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
#' df_scores<-find_scores(formula=Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   wlr="mw",
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
                       wlr,
                       t_star = NULL,
                       s_star = NULL,
                       rho = NULL,
                       gamma = NULL){
  
  w<-find_weights(formula=formula,
                 data=data,
                 wlr=wlr,
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
  #####################################
  ## include treatment group (DM)

  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  group_col<-formula_vars[-(1:2)]

  data$t_j<-data[[time_col]]
  data$event<-data[[status_col]]
  data$group<-data[[group_col]]

  df<-merge(x = df[,c("t_j","score_censored","score_event")], y = data[,c("t_j","event","group")], by = "t_j", all.x = TRUE)
  ###########################################
  ## pick out score (event or censored) (DM)
  df$score <- with(df, event * score_event + (1 - event) * score_censored)
  ## standardize scores to (-1, 1) (DM)
  max_a <- max(df$score)
  min_a <- min(df$score)
  A = 2 / (max_a - min_a)
  B = 1 - A * max_a
  df$standardized_score <- df$score * A + B
  ############################################
  
  df$rank<-factor(df$t_j)
  levels(df$rank)<-paste0("(",1:length(unique(df$t_j)),")")
  out<-list(df=df)
  class(out)<-"wlrt_score"
  out
}



#' @export
plot.wlrt_score<-function(x,...){
  df<-x$df
  
  ## suggested changes (DM)
  #df$x_pos<-1:length(unique(df$rank))
  df$x_pos<-1:length(df$rank)
  
  # df_experimental_cens<-df[df$group == "experimental" & df$event == 0,]
  # df_experimental_event<-df[df$group == "experimental" & df$event == 1,]
  # df_control_cens<-df[df$group == "control" & df$event == 0,]
  # df_control_event<-df[df$group == "control" & df$event == 1,]
  
  group_labels <- unique(df$group)
  gl1 <- group_labels[1]
  gl2 <- group_labels[2]
  
  df_gl1_cens<-df[df$group == gl1 & df$event == 0,]
  df_gl1_event<-df[df$group == gl1 & df$event == 1,]
  df_gl2_cens<-df[df$group == gl2 & df$event == 0,]
  df_gl2_event<-df[df$group == gl2 & df$event == 1,]
  
  args <- list(ylim=c(-1,1))
  
  inargs <- list(...)
  args[names(inargs)] <- inargs
  
  do.call(plot,c(list(x=df$x_pos,
                      type="n",xlab="Rank",xaxt="n",ylab="Score"),args))
  
  axis(side=1, labels=df$rank, at=df$x_pos)
  
  mycol_blue <- rgb(0, 0, 255, max = 255, alpha = 100)
  mycol_red <- rgb(255, 0, 0, max = 255, alpha = 100)
  #points(x=df_control_event$x_pos,y=df_control_event$standardized_score,pch=15,col=mycol_red,cex=1)
  #points(x=df_control_cens$x_pos,y=df_control_cens$standardized_score,pch=0,col=mycol_red,cex=1)
  #points(x=df_experimental_event$x_pos,y=df_experimental_event$standardized_score,pch=16,col=mycol_blue,cex=1)
  #points(x=df_experimental_cens$x_pos,y=df_experimental_cens$standardized_score,pch=1,col=mycol_blue,cex=1)
  points(x=df_gl1_event$x_pos,y=df_gl1_event$standardized_score,pch=15,col=mycol_red,cex=1)
  points(x=df_gl1_cens$x_pos,y=df_gl1_cens$standardized_score,pch=0,col=mycol_red,cex=1)
  points(x=df_gl2_event$x_pos,y=df_gl2_event$standardized_score,pch=16,col=mycol_blue,cex=1)
  points(x=df_gl2_cens$x_pos,y=df_gl2_cens$standardized_score,pch=1,col=mycol_blue,cex=1)
}
