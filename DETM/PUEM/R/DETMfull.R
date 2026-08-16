#' DETM fitting with full parameters
#'
#' @param posdata n by p feature matrix for the source (positive) data in PU data.
#' @param negdata m by p feature matrix for the target (unlabled) data in PU data.
#' @param pival   initial value of the mixture proportion of positive data in the target dataset, is set for the EM algorithm.
#' @param alp1     initial value of alpha_1 of DETM in EM algorithm
#' @param beta1    A vector of initial value of beta_1 of DETM in EM algorithm
#' @param alp2     initial value of alpha_2 of DETM in EM algorithm
#' @param beta2    A vector of initial value of beta_2 of DETM in EM algorithm
#' @param maxite     maximum iteration number of EM algorithm
#' @param eps     tolerance value about stopping the algorithm when the increment in the log-EL after an iteration is no greater than, e.g. 1e-4.
#' @return  A list including mixture proportion estimation, parameter estimation, iterations, likelihood and penalized likelihood in DETM.
#' @export
#'
DETMfull=function(posdata,negdata,pival,alp1,beta1,alp2,beta2,maxite=1000,eps=1e-4)
{

  n=nrow(posdata)
  m=nrow(negdata)
  d=ncol(negdata)
  Lambda<-1/(n+m)^(3/2)
  itenum=0
  err=1
  newlik<- -10^30
  newlik_penal<- -10^30
  itenum<-0
  if(pival<1/m){pival=1/m}
  if(pival>1-1/m){pival=1-1/m}
  while(err>eps & itenum<=maxite)
  {
    old_itenum=itenum
    itenum<-itenum+1
    old_alp1<-alp1
    old_alp2<-alp2
    old_beta1<-beta1
    old_beta2<-beta2
    old_pival<-pival
    oldlik_penal=newlik_penal
    old_newlik<-newlik
    ti=old_alp2-old_alp1+negdata%*%(old_beta2-old_beta1)
    ti=as.numeric(ti)

    wi=old_pival/(old_pival+(1-old_pival)*exp(ti))

    pival=mean(wi)
    if(pival<1/m){pival=1/m}
    if(pival>1-1/m){pival=1-1/m}

    dlab=c(rep(0,n),rep(1,m),rep(2,m))

    weight=c(rep(1,n),wi,1-wi)

    xmat=rbind(posdata,negdata,negdata)
    ind_break <- FALSE
    ind_break <- tryCatch(
      expr = {
        out=glmnet::glmnet(x=xmat,y=dlab,family="multinomial",weights=weight,alpha=0,lambda=Lambda,standardize = FALSE)
      },
      error = function(e){
        print("ERROR, from DETMfull")
        return(TRUE)
      },
      warning = function(w){
        print(" Warning, from DETMfull")
        return(TRUE)
      })

    if(isTRUE(ind_break)){
      return(list( penal_loglik=NA))
    }

    out=coef(out)

    theta=rbind(as.numeric(out$`1`-out$`0`), as.numeric(out$`2`-out$`0`) )

    salp1=theta[1,1]
    beta1=as.numeric(theta[1,-1])
    salp2=theta[2,1]
    beta2=as.numeric(theta[2,-1])

    hatpi1=1/(1+exp(salp1+posdata%*%beta1)+exp(salp2+posdata%*%beta2))
    hatpi2=1/(1+exp(salp1+negdata%*%beta1)+exp(salp2+negdata%*%beta2))

    hatpi=c(as.numeric(hatpi1),as.numeric(hatpi2))
    n_hatpi=hatpi/sum(hatpi)

    alp1=salp1-log(sum(wi)/n)
    alp2=salp2-log(sum(1-wi)/n)


    part1=sum(log(n_hatpi[1:n]))
    part2<- sum(as.numeric(log(1/(1+1/(exp(salp1+negdata%*%beta1)
                                       +exp(salp2+negdata%*%beta2))))))
    newlik<-part1+part2+m*log(n/(sum(hatpi)*m))

    obeta0<-as.numeric(out$`0`)
    obeta1<-as.numeric(out$`1`)
    obeta2<-as.numeric(out$`2`)
    newlik_penal=newlik-(1/2*Lambda)*(obeta0[-1]%*%obeta0[-1]+obeta1[-1]%*%obeta1[-1]+obeta2[-1]%*%obeta2[-1])
    err=newlik_penal-oldlik_penal

  }
  if(err<0){
    return(
      list(
        pival=old_pival,
        alp1=old_alp1,
        beta1=old_beta1,
        alp2=old_alp2,
        beta2=old_beta2,
        loglik=old_newlik,itenum=old_itenum,
        penal_loglik=oldlik_penal
      )
    )

  }
  else{
    return(
      list(
        pival=pival,
        alp1=alp1,
        beta1=beta1,
        alp2=alp2,
        beta2=beta2,
        loglik=newlik,itenum=itenum,
        penal_loglik=newlik_penal
      )
    )}
}
