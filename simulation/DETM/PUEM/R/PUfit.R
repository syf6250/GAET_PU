#' PU fitting with four models (SETMnull;SETMfull;DETMnull;DETMfull)
#'
#' @param posdata n by p feature matrix for the source (positive) data in PU data.
#' @param negdata m by p feature matrix for the target (unlabled) data in PU data.
#' @param pival   initial value of the mixture proportion of positive data in the target dataset, is set for the EM algorithm.
#' @param alp1     initial value of alpha_1 in EM algorithm
#' @param beta1    A vector of initial value of beta_1 in EM algorithm
#' @param alp2     initial value of alpha_2 in EM algorithm
#' @param beta2    A vector of initial value of beta_2 in EM algorithm
#' @return  A list including four models with their mixture proportion estimation, parameter estimation, iterations, likelihood and penalized likelihood.
#' @export
#'
PUfit=function(posdata,negdata,pival,alp1,beta1,alp2,beta2)
{

  out_detm_full<-list(
    penal_loglik=-10^30
  )
  out_setm_full<-out_detm_full
  out_detm_null<-out_detm_full
  out_setm_null<-out_detm_full

  d<-ncol(negdata)

  out2_full<-SETMfull(posdata,negdata,pival,alp2,beta2)
  out2_null<-SETMnull(posdata,negdata,pival,alp2,beta2)
  out3_full<-DETMfull(posdata,negdata,pival,alp1,beta1,alp2,beta2)
  out3_null<-DETMnull(posdata,negdata,pival,alp1,beta1,alp2,beta2)

  if(!is.na(out2_full$penal_loglik)){
    out_setm_full<-out2_full }
  if(!is.na(out2_null$penal_loglik)){
    out_setm_null<-out2_null   }
  if(!is.na(out3_full$penal_loglik)){
    out_detm_full<-out3_full}
  if(!is.na(out3_null$penal_loglik)){
    out_detm_null<-out3_null}


  alp1_try<-alp1+runif(9,-0.5,0.5)
  beta1_try<-beta1+matrix(runif(9*d,-0.5,0.5),nrow=9)
  alp2_try<-alp2+runif(9,-0.5,0.5)
  beta2_try<-beta2+matrix(runif(9*d,-0.5,0.5),nrow=9)
  if(pival>0.05 & pival<0.95) {pival_try<-pival+runif(9,-0.15,0.15)}
  if(pival<=0.05){pival_try<-pival+runif(9,0.035,0.1)}
  if(pival>=0.95){pival_try<-pival+runif(9,-0.1,-0.035)}


  for(j in 1:9){
    out2_new_null<-SETMnull(posdata,negdata,pival,alp2_try[j],beta2_try[j,])
    if(is.na(out2_new_null$penal_loglik)){next}
    if(out2_new_null$penal_loglik>=out_setm_null$penal_loglik){out_setm_null<-out2_new_null}
  }

  for (j in 1:8){
    out2_new_full<-SETMfull(posdata,negdata,pival_try[j],alp2_try[j],beta2_try[j,])
    if(is.na(out2_new_full$penal_loglik)){next}
    if(out2_new_full$penal_loglik>=out_setm_full$penal_loglik){out_setm_full<-out2_new_full} # record the largest likelihood
  }
  if(!is.na(out_setm_null$penal_loglik)){
    out2_new_full<-SETMfull(posdata,negdata,pival,out_setm_null$alp,out_setm_null$beta)
    if(!is.na(out2_new_full$penal_loglik)){
      if(out2_new_full$penal_loglik>=out_setm_full$penal_loglik){out_setm_full<-out2_new_full }}}

  for (j in 1:8){
    out3_new_null<-DETMnull(posdata,negdata,pival,alp1_try[j],beta1_try[j,],alp2_try[j],beta2_try[j,])
    if(is.na(out3_new_null$penal_loglik)){next}
    if(out3_new_null$penal_loglik>=out_detm_null$penal_loglik){out_detm_null<-out3_new_null }
  }
  if(!is.na(out_setm_full$penal_loglik)){
    out3_new_null<-DETMnull(posdata,negdata,pival,0,rep(0,d),out_setm_full$alp,out_setm_full$beta)
    if(!is.na(out3_new_null$penal_loglik)){
      if(out3_new_null$penal_loglik>=out_detm_null$penal_loglik){out_detm_null<-out3_new_null }}}

  for (j in 1:7){
    out3_new_full<-DETMfull(posdata,negdata,pival_try[j],alp1_try[j],beta1_try[j,],alp2_try[j],beta2_try[j,])
    if(out3_new_full$penal_loglik>=out_detm_full$penal_loglik){out_detm_full<-out3_new_full }
  }

  if(!is.na(out_detm_null$penal_loglik)){
    out3_new_full<-DETMfull(posdata,negdata,out_detm_null$pival,out_detm_null$alp1,out_detm_null$beta1,out_detm_null$alp2,out_detm_null$beta2)
    if(!is.na(out3_new_full$penal_loglik)){
      if(out3_new_full$penal_loglik>=out_detm_full$penal_loglik){out_detm_full<-out3_new_full }}}

  if(!is.na(out_setm_full$penal_loglik)){
    out3_new_full<-DETMfull(posdata,negdata,out_setm_full$pival,0,rep(0,d),out_setm_full$alp,out_setm_full$beta)
    if(!is.na(out3_new_full$penal_loglik)){
      if(out3_new_full$penal_loglik>=out_detm_full$penal_loglik){out_detm_full<-out3_new_full}}}


  list(out_setm_null=out_setm_null,
       out_setm_full=out_setm_full,
       out_detm_null=out_detm_null,
       out_detm_full=out_detm_full)

}
