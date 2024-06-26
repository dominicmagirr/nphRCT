#' Calculate at-risk table
#' 
#' This function calculates the number of individuals at risk and number of events at each distinct event time
#' (and censoring time if `include_cens==TRUE`). 
#' 
#' @template formula
#' @template data
#' @template include_cens
#' @param timefix Deal with floating point issues (as in the survival package). Default is TRUE. May need to set FALSE for simulated data.
#' @return Data frame 
#' 
#' @return
#' Returns a data frame with the following columns:
#' - time `t_j` 
#' - number of events in each of the treatments at `t_j`
#' - combined number of events in both treatments at event time `t_j`
#' - number of individuals at risk in each of the treatment groups just before time `t_j`
#' - combined number of individuals at risk in both treatment groups just before time `t_j`
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
#'   n_c=5,
#'   n_e=5,
#'   max_cal_t = 36
#' )
#' #with censoring times included
#' find_at_risk(formula=survival::Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   include_cens=TRUE)
#' #with censoring times excluded
#' find_at_risk(formula=survival::Surv(event_time,event_status)~group,
#'   data=sim_data,
#'   include_cens=FALSE)
#' @export

find_at_risk<-function(formula,
                       data,
                       include_cens=TRUE,
                       timefix = TRUE){
  
  
  
  
  check_formula(formula=formula,
                data=data)
  
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  terms_vars<-formula_vars[-(1:2)]
  Terms <- stats::terms(formula,"strata")
  strata_index <- survival::untangle.specials(Terms,"strata")$terms
  if (length(strata_index)>0){stop("Function does not account for strata")}

  #Find at-risk tables
  group_col<-terms_vars
  groups <- as.vector(sort(unique(data[[group_col]])))
  if (length(groups)!=2){stop("Only 2 treatment groups allowed")}
  
  data<-data[order(data[[time_col]]),]
  trt <- data[[group_col]]
  times <- data[[time_col]]
  d_j <- data[[status_col]]
  
  
  ################
  ## timefix logic
  ################
  if (timefix){
  times <- survival::aeqSurv(survival::Surv(times, rep(1, length(times))))[,1]
  }
  
  if (include_cens==TRUE){
    #number of events
    n_event_g <- matrix(0, length(times), 2)
    for (i in 1:2){
      n_event_g[, i] <- (trt == groups[i])*d_j
    }
    n_event_g<-data.frame(times,n_event_g)
    n_event_g<-stats::aggregate(. ~ times, n_event_g, FUN = sum)[,-1]
    
    #number at risk
    t_j<-unique(times)
    n_risk_g <- matrix(0, length(t_j), 2)
    for (i in 1:2){
      n_risk_g[, i] <- colSums(matrix(rep(t_j,each = sum(trt == groups[i])),ncol=length(t_j)) <= times[trt == groups[i]])
    }
  }else{
    #number of events
    d_j_events<-(d_j > 0)
    trt_events <- trt[d_j_events]
    t_j_events<-unique(times[d_j_events])
    n_event_g <- matrix(table(trt[d_j_events > 0], times[d_j_events]),ncol=2,byrow = TRUE)
    #number at risk
    n_risk_g <- matrix(0, length(t_j_events), 2)
    for (i in 1:2){
      n_risk_g[, i] <- 
        colSums(matrix(rep(t_j_events,each = sum(trt == groups[i])),ncol=length(t_j_events)) <= times[trt == groups[i]])
    }
  }
  n_event <- rowSums(n_event_g)
  n_risk <- rowSums(n_risk_g)

  out<-data.frame(n_event_g,n_event,n_risk_g,n_risk)
  colnames(out)<-c(paste0("n_event_",groups),"n_event",paste0("n_risk_",groups),"n_risk")
  
  if (include_cens==TRUE){
    out<-cbind(t_j=t_j,out)
  }else{
    out<-cbind(t_j=t_j_events,out)
  }
  
  out

}
