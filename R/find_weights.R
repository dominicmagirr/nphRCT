#' @export

find_weights<-function(data,
                           formula,
                           wlr,
                           t_star = NULL,
                           s_star = NULL,
                           rho = NULL,
                           gamma = NULL,
                       include_cens=TRUE){
  formula_vars <- all.vars(formula)
  time_col <- formula_vars[1]
  status_col <- formula_vars[2]

  #Find weights
  wlr <- match.arg(wlr, c("lr", "fh", "mw"))

  Surv <- survival::Surv
  Y <- eval(formula[[2]], data) ##formula[[2]] takes lhs

  km_fit <- survival::survfit(Surv(get(time_col),get(status_col))~1,
                                data = data,
                                timefix = FALSE)
  if(include_cens==TRUE){
    t_j <- km_fit$time
    S_hat <- km_fit$surv
    S_hat_minus <- c(1, S_hat[1:(length(S_hat) - 1)])
  }else{
    t_j <- km_fit$time[km_fit$n.event > 0]
    S_hat <- km_fit$surv[km_fit$n.event > 0]
    S_hat_minus <- c(1, S_hat[1:(length(S_hat) - 1)])
  }

  if (wlr == "lr") {
    w <- rep(1, length(S_hat_minus))
  }
  if (wlr == "fh") {
    if (is.null(rho) ||
        is.null(gamma))
      stop("Must specify rho and gamma")
    if (rho < 0 ||
        gamma < 0)
      stop("rho and gamma must be non-negative")

    w <- S_hat_minus ^ rho * (1 - S_hat_minus) ^ gamma
  }
  if (wlr == "mw") {
    if (is.null(t_star) && is.null(s_star)) stop("must specify either t_star or s_star")
    if (!is.null(t_star) && !is.null(s_star)) stop("must specify either t_star or s_star (not both)")

    if (!is.null(t_star)){
      if (any(t_j < t_star)) {
        w <- pmin(1 / S_hat_minus,
                  1 / S_hat[max(which(t_j < t_star))])
      }
      else {
        w <- rep(1, length(S_hat_minus))
      }
    }
    else {
      w <- pmin(1 / S_hat_minus,
                1 / s_star)
    }
  }

  out<-list(w=w,km_fit=km_fit)
  out
}
