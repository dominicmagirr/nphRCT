find_pseudovalues <- function(formula,
                              data,
                              method = "RMST",
                              tau = 1){
  
  ##################################################
  ### rmst1 code from survRM2 package
  rmst1=function(time, status, tau, alpha=0.05){
    #-- time
    #-- statuts
    #-- tau -- truncation time
    #-- alpha -- gives (1-alpha) confidence interval
    
    ft= survfit(Surv(time, status)~1)
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
    out[1,]=c(rmst, rmst.se, rmst-qnorm(1-alpha/2)*rmst.se, rmst+qnorm(1-alpha/2)*rmst.se)
    out[2,]=c(tau-out[1,1], rmst.se, tau-out[1,4], tau-out[1,3])
    rownames(out)=c("RMST","RMTL")
    colnames(out)=c("Est.", "se", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
    
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
  #############################################################
  ### extract terms from formula
  check_formula(formula=formula,data=data)
  Surv<-survival::Surv
  
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]
  terms_vars<-formula_vars[-(1:2)]
  Terms <- terms(formula,"strata")
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
  
  
  rmst_full <- rmst1(time = data[[time_col]],
                     status = data[[status_col]],
                     tau = tau)$rmst[1]
  
  n <- length(data[[time_col]])
  
  ### also need the size of each arm...
  
  p_v_rmst <- numeric(length(dat[[time_col]]))
  
  for (i in seq_along(p_v_rmst)){
    
    data_minus_i <- data[-i,]
    
    rmst_minus_i <- rmst1(time = data_minus_i[[time_col]],
                          status = data_minus_i[[status_col]],
                          tau = tau)$rmst[1]
    
    p_v_rmst[i] <- n * rmst_full - (n - 1) * rmst_minus_i
    
  }
  
  max_p <- max(-p_v_rmst)
  min_p <- min(-p_v_rmst)
  
  A = 2 / (max_p - min_p)
  B = 1 - A * max_p
  
  df_rmst_pseudo <- data.frame(time = data[[time_col]],
                               event = data[[status_col]],
                               arm = data[[group_col]],
                               scores = -p_v_rmst * A + B)
  
  
  df_rmst_pseudo
  
  
  
}