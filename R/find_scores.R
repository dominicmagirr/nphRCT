#' @export

find_scores<-function(data,
                           formula,
                           wlr,
                           t_star = NULL,
                           s_star = NULL,
                           rho = NULL,
                           gamma = NULL){

  out_weights<-find_weights(data=data,
                 formula=formula,
                 wlr=wlr,
                 t_star = t_star,
                 s_star = s_star,
                 rho = rho,
                 gamma = gamma,
                 include_cens=TRUE)
  w<-out_weights$w
  km_fit<-out_weights$km_fit

  df<-data.frame(time=km_fit$time,
                 n_risk=km_fit$n.risk,
                 n_event=km_fit$n.event,
                 n_censor=km_fit$n.censor,
                 w=w,
                 S_hat=km_fit$surv
  )
  df$S_hat_minus<-with(df,c(1, S_hat[1:(length(S_hat) - 1)]))

  df$score_cens<-with(df,-cumsum(w*(n_event/n_risk)))
  df$score_event<-with(df,score_cens+w)
  df$rank<-paste0("(",nrow(df):1,")")
  out<-list(scores_df=df)
  class(out) <- "wlr_scores"
  out

}

#' @export
plot.wlr_scores<-function(x,...){
  scores_df<-x$scores_df
  scores_df$x_pos<-1:nrow(scores_df)
  scores_df_cens<-scores_df[scores_df$n_censor>0,]
  scores_df_event<-scores_df[scores_df$n_event>0,]

  args <- list(ylim=c(min(scores_df_cens$score_cens,scores_df_event$score_event)-0.5,
                       max(scores_df_cens$score_cens,scores_df_event$score_event)+0.5))
  inargs <- list(...)
  args[names(inargs)] <- inargs


  do.call(plot,c(list(x=scores_df$x_pos,
       type="n",xlab="Rank",xaxt="n",ylab="Score"),args))

  axis(side=1, labels=scores_df$rank, at=scores_df$x_pos)
  points(x=scores_df_cens$x_pos,y=scores_df_cens$score_cens,pch=21,bg="white")
  points(x=scores_df_event$x_pos,y=scores_df_event$score_event,pch=21,bg="black")
}
