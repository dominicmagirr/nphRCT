#' Weighted log-rank test
#'
#' This function can perform two types of weighted log-rank test,
#' the modestly-weighted log-rank test and the Fleming-Harrington (\eqn{\rho},\eqn{\gamma}) test, in addition to the standard log-rank test.
#'
#' @template data
#' @template formula
#' @template wlr
#' @template t_star
#' @template s_star
#' @template rho
#' @template gamma
#' @return List containing the outcome of the weighted log-rank test.
#' `u` is the test statistic U for the weighted log-rank test
#' `v_u` is the variance of test statistic U
#' `z` is the Z-score
#' `trt_group` indicates which of the treatment arms the test statistic U corresponds to
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
#' #example setting t_star
#' wlrt(data=sim_data,
#'   formula=Surv(event_time,event_status)~group,
#'   wlr="mw",
#'   t_star = 3
#' )
#' #example setting s_star
#' wlrt(data=sim_data,
#'   formula=Surv(event_time,event_status)~group,
#'   wlr="mw",
#'   s_star = 0.5
#' )
#'
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

wlrt <- function(data,
                 formula,
                 wlr,
                 t_star = NULL,
                 s_star = NULL,
                 rho = NULL,
                 gamma = NULL) {

  #Checks
  Surv<-survival::Surv
  Y <- eval(formula[[2]], data) ##formula[[2]] takes lhs
  if (!inherits(Y, "Surv")){
    stop("formula has an incorrect format")
  }
  if (!attr(Y, "type") == "right") {
    stop("Censoring type should be right censoring")
  }

  Terms <- terms(formula,"strata")
  terms_vars <- labels(Terms)
  strata_vars <- attr(Terms,"specials")$strata

  #no strata:
  if (is.null(strata_vars)) {
    if (length(terms_vars)!=1) {
      stop("One treatment arm indicator must be specified in formula")
    }
    if (any(is.na(data[, formula_vars]))) {
      stop("NA's in data set. wlrt doesn't have a default for missing data.")
    }

    #calculate weights
    w<-find_weights(data=data,
                 formula=formula,
                 wlr=wlr,
                 t_star = t_star,
                 s_star = s_star,
                 rho = rho,
                 gamma = gamma,
                 include_cens=FALSE)$w

    #at risk tables
    df_events<-find_at_risk(data=data,
                            formula=formula,
                           include_cens=FALSE)
    n_eventg <- as.matrix(df_events[, 4:5])
    n_event<-df_events$n_event
    n_riskg <- as.matrix(df_events[, 2:3])
    n_risk<-df_events$n_risk

    # test statistics
    observed <- w %*% n_eventg
    expected <- w %*% (n_riskg * (n_event / n_risk))
    u <- (observed - expected)[, 2]

    v_u <-w ^ 2 * (apply(n_riskg, 1, prod) * n_event * (n_risk - n_event)) / (n_risk ^ 2 * (n_risk - 1))
    v_u[n_risk == 1] <- 0
    v_u<-sum(v_u)

    z = u / sqrt(v_u)

    groups <- sort(unique(df_events$group))

    out <- list(
      u = as.vector((observed - expected)[, 2]),
      v_u = as.vector(v_u),
      z = as.vector(z),
      trt_group = sort(unique(data[[terms_vars]]))[2]
    )
    out
  }

  #strata
  else{

  }

}

