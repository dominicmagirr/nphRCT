#' Calculate scores
#'
#' Weighted log-rank tests can also be thought in terms of assigning a score to the each of
#' the events (including censoring) and comparing the average score on each arm, see Magirr (2021) <doi:10.1002/pst.2091>. 
#' This function calculates the scores
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
#' @param timefix Deal with floating point issues (as in the survival package). Default is TRUE. May need to set FALSE for simulated data.
#' @return 
#' Data frame. Each row corresponds to an event or censoring time.
#' At each time specified in `t_j` the columns indicate
#' - `event` the event indicator
#' - `group` the treatment arm indicator
#' - `score` the assigned score at time `t_j`
#' - `standardized_score` the value of `score` standardized to be between -1 and 1
#' 
#' @details
#'
#' Select which of the tests to perform using argument `method`.
#' For the weighted log-rank tests, the output is calculated as outlined in `vignette("weighted_log_rank_tests", package="nphRCT")`.
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
#' <doi:10.1002/pst.2091>
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
                      tau=NULL,
                      timefix = TRUE){
  
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
                   include_cens=TRUE,
                   timefix = timefix)
    
    df <- find_at_risk(formula=formula,
                       data=data,
                       include_cens=TRUE,
                       timefix = timefix)
    
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
  
  df[["event"]]<-factor(df[["event"]],levels=c(1,0))
  levels(df[["event"]])<-c( "event","censored")
  
  df[["event_group"]]<-factor(paste(df[["group"]],df[["event"]],sep=", "),
                              levels=c(paste(gl1,"event",sep=", "),paste(gl2,"event",sep=", "),
                                       paste(gl1,"censored",sep=", "),paste(gl2,"censored",sep=", ")))
  
  args <- list(col="",x="Time",y="Score")
  inargs <- list(...)
  args[names(inargs)] <- inargs
  labels<-do.call(ggplot2::labs,args)
  
  ggplot2::ggplot(df, ggplot2::aes_string(x="t_j", y="standardized_score",col="event_group")) + ggplot2::geom_point() +
    ggplot2::ylim(-1.1,1.1)+
    labels+ ggplot2::scale_color_manual(values = c("#F8766D", "#00BFC4", "lightsalmon", "darkslategray2"))+
    ggplot2::geom_hline(yintercept = mean(df[df[["group"]]==gl1,"standardized_score"]), color="#F8766D",linetype=2)+
    ggplot2::geom_hline(yintercept = mean(df[df[["group"]]==gl2,"standardized_score"]), color="#00BFC4",linetype=2)+
    ggplot2::theme_classic()
  }
