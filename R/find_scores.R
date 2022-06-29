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
#' @return Data frame containing output from function `find_weights` and also additional columns:
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
  
  #####################################
  ## include treatment group (DM)
  df_at_risk <- find_at_risk(formula=formula,
                             data=data,
                             include_cens=TRUE)
  
  df_at_risk$n_censored_control <- rev(diff(rev(c(df_at_risk$n_risk_control,0))) - rev(df_at_risk$n_event_control))
  df_at_risk$n_censored_experimental <- rev(diff(rev(c(df_at_risk$n_risk_experimental,0))) - rev(df_at_risk$n_event_experimental))
  
  df$n_event_control <- df_at_risk$n_event_control
  df$n_event_experimental <- df_at_risk$n_event_experimental
  df$n_censored_control <- df_at_risk$n_censored_control
  df$n_censored_experimental <- df_at_risk$n_censored_experimental
  
  ###########################################
  ## replicate rows with ties (DM)
  df_full <- NULL
  for (i in seq_along(df$t_j)) {
    
    df_i_e_e <- purrr::map_df(rep(i, df[i,"n_event_experimental"]), function(x) df[x,])
    df_i_e_e$group <- "experimental"
    df_i_e_e$event <- 1
    
    df_i_c_e <- purrr::map_df(rep(i, df[i,"n_event_control"]), function(x) df[x,])
    df_i_c_e$group <- "control"
    df_i_c_e$event <- 1
    
    df_i_e_c <- purrr::map_df(rep(i, df[i,"n_censored_experimental"]), function(x) df[x,])
    df_i_e_c$group <- "experimental"
    df_i_e_c$event <- 0
    
    df_i_c_c <- purrr::map_df(rep(i, df[i,"n_censored_control"]), function(x) df[x,])
    df_i_c_c$group <- "control"
    df_i_c_c$event <- 0
    
    df_full <- rbind(df_full, df_i_e_e, df_i_c_e, df_i_e_c, df_i_c_c)
  }
  df <- df_full
  ###########################################
  ## pick out score (event or censored) (DM)
  df$score <- with(df, event * score_event + (1 - event) * score_cens)
  ## standardize scores to (-1, 1) (DM)
  max_a <- max(df$score)
  min_a <- min(df$score)
  A = 2 / (max_a - min_a)
  B = 1 - A * max_a
  df$standardized_score <- df$score * A + B
  ############################################
  
  
  
  df$rank<-paste0("(",nrow(df):1,")")

  out<-list(df=df)
  class(out)<-"wlrt_score"
  out
}



#' @export
plot.wlrt_score<-function(x,...){
  df<-x$df
  df$x_pos<-1:nrow(df)
  #df_cens<-df[df$n_censor>0,]
  #df_event<-df[df$n_event>0,]
  df_experimental_cens<-df[df$group == "experimental" & df$event == 0,]
  df_experimental_event<-df[df$group == "experimental" & df$event == 1,]
  df_control_cens<-df[df$group == "control" & df$event == 0,]
  df_control_event<-df[df$group == "control" & df$event == 1,]
  
  #args <- list(ylim=c(min(df_cens$score_cens, df_event$score_event)-0.5,
  #                    max(df_cens$score_cens, df_event$score_event)+0.5))
  
  args <- list(ylim=c(-1,1))
  
  inargs <- list(...)
  args[names(inargs)] <- inargs
  
  do.call(plot,c(list(x=df$x_pos,
                      type="n",xlab="Rank",xaxt="n",ylab="Score"),args))
  
  axis(side=1, labels=df$rank, at=df$x_pos)
  
  #points(x=df_cens$x_pos,y=df_cens$score_cens,pch=21,bg="white")
  #points(x=df_event$x_pos,y=df_event$score_event,pch=21,bg="black")
  
  points(x=df_experimental_cens$x_pos,y=df_experimental_cens$standardized_score,pch=21,bg="white")
  points(x=df_experimental_event$x_pos,y=df_experimental_event$standardized_score,pch=21,bg="black")
  points(x=df_control_cens$x_pos,y=df_control_cens$standardized_score,pch=21,bg="white", col = 2)
  points(x=df_control_event$x_pos,y=df_control_event$standardized_score,pch=21,bg="red", col = 2)
  
}
