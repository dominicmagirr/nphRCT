#' @export

find_at_risk<-function(data,
                       formula,
                       include_cens=TRUE){

  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  Terms <- terms(formula,"strata")
  terms_vars <- labels(Terms)
  strata_vars <- attr(Terms,"specials")$strata
  if (length(strata_vars)>0){stop("There should be no strata in formula argument for function find_at_risk")}

  #Find at-risk tables
  group_col<-terms_vars
  df_events<-data.frame(time=data[[time_col]],
                        status=data[[status_col]],
                        group=data[[group_col]]
  )
  groups <- sort(unique(df_events$group))
  if (length(groups)!=2){stop("Only 2 treatment groups allowed")}
  df_events[,paste0("n_risk_", groups)]<-sapply(groups,function(i){+(df_events$group == i)})
  df_events[,paste0("n_event_", groups)]<-sapply(groups,function(i){df_events$status * df_events[[paste0("n_risk_", i)]]})
  df_events[["group"]] <- df_events[["status"]] <-NULL
  df_events <- aggregate(. ~ time, df_events, FUN = sum)
  df_events[,paste0("n_risk_", groups)]<-sapply(groups,function(i){rev(cumsum(rev(df_events[,paste0("n_risk_", i)])))})

  if(include_cens==FALSE){
    n_event<-rowSums(df_events[,paste0("n_event_", groups)])
    df_events<-df_events[n_event>0,]
  }

  n_eventg <- as.matrix(df_events[, paste0("n_event_", groups)])
  df_events[["n_event"]] <- rowSums(n_eventg)
  n_riskg <- as.matrix(df_events[, paste0("n_risk_", groups)])
  df_events[["n_risk"]] <- rowSums(n_riskg)

  df_events

}
