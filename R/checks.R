check_formula<-function(formula,data){
  Surv<-survival::Surv
  Y <- eval(formula[[2]], data) ##formula[[2]] takes lhs
  if (!inherits(Y, "Surv")){
    stop("formula has an incorrect format")
  }
  if (!attr(Y, "type") == "right") {
    stop("Censoring type should be right censoring")
  }
  formula_vars <- all.vars(formula)
  if (any(is.na(data[, formula_vars]))) {
    stop("NAs in data set, no default for missing data.")
  }  
}

check_lr<-function(t_star,s_star,rho,gamma){
  if (!is.null(t_star) || !is.null(s_star)) stop("do not specify t_star or s_star for log rank test 'lr'")
  if (!is.null(rho) || !is.null(gamma)) stop("do not specify rho or gamma for log rank test 'lr'")
}

check_fh<-function(t_star,s_star,rho,gamma){
  if (is.null(rho) ||
      is.null(gamma))
    stop("Must specify rho and gamma")
  if (rho < 0 ||
      gamma < 0)
    stop("rho and gamma must be non-negative")
  if (!is.null(t_star) || !is.null(s_star)) stop("do not specify t_star or s_star for Fleming-Harrington test 'fh'")
  
}

check_mw<-function(t_star,s_star,rho,gamma){
  if (is.null(t_star) && is.null(s_star)) stop("must specify either t_star or s_star")
  if (!is.null(t_star) && !is.null(s_star)) stop("must specify either t_star or s_star")
  if (!is.null(rho) || !is.null(gamma)) stop("do not specify rho or gamma for modestly weighted log rank test 'mw'")
}

