### rmst1 code from survRM2 package
rmst1<-function(time, status, tau, alpha=0.05){
  #-- time
  #-- statuts
  #-- tau -- truncation time
  #-- alpha -- gives (1-alpha) confidence interval
  Surv<-survival::Surv
  
  ft= survival::survfit(Surv(time, status)~1)
  idx=ft$time<=tau
  
  wk.time=sort(c(ft$time[idx],tau))
  wk.surv=ft$surv[idx]
  wk.n.risk =ft$n.risk[idx]
  wk.n.event=ft$n.event[idx]
  
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst = sum(areas)
  rmst
  
  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var =c(wk.var,0)
  rmst.var = sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  = sqrt(rmst.var)
  
  #--- check ---
  # print(ft, rmean=tau)
  
  #--- output ---
  out=matrix(0,2,4)
  out[1,]=c(rmst, rmst.se, rmst-stats::qnorm(1-alpha/2)*rmst.se, rmst+stats::qnorm(1-alpha/2)*rmst.se)
  out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
  rownames(out)=c("RMST","RMTL")
  colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), 
                  paste("upper .",round((1-alpha)*100, digits=0), sep=""))
  
  Z=list()
  Z$result=out
  Z$rmst = out[1,]
  Z$rmtl = out[2,]
  Z$tau=tau
  Z$rmst.var = rmst.var
  Z$fit=ft
  class(Z)="rmst1"
  
  return(Z)
  
}

find_pseudovalues <- function(formula,
                              data,
                              method,
                              tau = 1){

  method<-match.arg(method,c("rmst","ms"))
  
  ### extract terms from formula
  check_formula(formula=formula,data=data)
  
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
  Surv<-survival::Surv
  
  if(method=="rmst"){
    
    rmst_full <- rmst1(time = data[[time_col]],
                       status = data[[status_col]],
                       tau = tau)$rmst[1]
    n <- nrow(data)
    
    ### also need the size of each arm...
    
    p_v_rmst <- numeric(nrow(data))
    
    for (i in seq_along(p_v_rmst)){
      data_minus_i <- data[-i,]
      
      rmst_minus_i <- rmst1(time = data_minus_i[[time_col]],
                            status = data_minus_i[[status_col]],
                            tau = tau)$rmst[1]
      
      p_v_rmst[i] <- n * rmst_full - (n - 1) * rmst_minus_i
      
    }
  
    df_rmst_pseudo <- data.frame(t_j = data[[time_col]],
                                 event = data[[status_col]],
                                 group = data[[group_col]],
                                 score = -p_v_rmst )
    df_rmst_pseudo<-df_rmst_pseudo[order(df_rmst_pseudo[["t_j"]]),]

    return(df_rmst_pseudo)

    
  }
  if(method=="ms"){
    n <- length(data[[time_col]])
    
    ### non-parametric estimate of survival
    formula_km<-stats::as.formula(paste0("Surv(",time_col,",",status_col,")~1"))
    surv <- survival::survfit(formula_km, data = data)
    
    s_full <- summary(surv, time = tau)$surv
    
    p_v_surv <- numeric(length(data[[time_col]]))
    
    for (i in seq_along(p_v_surv)){
      
      df_minus_i <- data[-i,]
      
      ### non-parametric estimate of survival probability at tau
      formula_km<-stats::as.formula(paste0("Surv(",time_col,",",status_col,")~1"))
      surv_minus_i <- survival::survfit(formula_km, data = df_minus_i)
      s_full_minus_i <- summary(surv_minus_i, time = tau)$surv
      
      p_v_surv[i] <- n * s_full - (n - 1) * s_full_minus_i
      
    }
    df_milestone_pseudo <- data.frame(t_j = data[[time_col]],
                                 event = data[[status_col]],
                                 group = data[[group_col]],
                                 score = -p_v_surv )
    df_milestone_pseudo<-df_milestone_pseudo[order(df_milestone_pseudo[["t_j"]]),]
 
    return(df_milestone_pseudo)
 
    
  }
}

