#' SETM fitting with parameters given mixture proportion of positive data in the target dataset
#'
#' @param posdata n by p feature matrix for the source (positive) data in PU data.
#' @param negdata m by p feature matrix for the target (unlabled) data in PU data.
#' @param pival   fixed mixture proportion of positive data in the target dataset, is set for the EM algorithm.
#' @param alp     initial value of alpha of SETM in EM algorithm
#' @param beta    A vector of initial value of beta of SETM in EM algorithm
#' @param maxite     maximum iteration number of EM algorithm
#' @param eps     tolerance value about stopping the algorithm when the increment in the log-EL after an iteration is no greater than, e.g. 1e-4.
#' @return  A list including fixed mixture proportion estimation, parameter estimation, iterations, log-likelihood and penalized log-likelihood in SETM.
#' @export
#'
SETMnull=function(posdata,negdata,pival,alp,beta,maxite=1000,eps=1e-4)
{

  n=nrow(posdata)
  m=nrow(negdata)
  Lambda<-1/(n+m)^(3/2)
  itenum=0
  err=1
  newlik<- -10^30
  newlik_penal<--10^30

  while(err>eps & itenum<=maxite)
  {
    old_itenum=itenum
    itenum<-itenum+1
    old_alp<-alp
    old_beta<-beta
    old_pival<-pival
    oldlik_penal=newlik_penal
    old_newlik<-newlik

    ti=old_alp+negdata%*%old_beta
    ti=as.numeric(ti)
    wi=old_pival/(old_pival+(1-old_pival)*exp(ti))

    dlab=c(rep(0,n),rep(0,m),rep(1,m))

    weight=c(rep(1,n),wi,1-wi)

    xmat=rbind(posdata,negdata,negdata)

    ind_break <- FALSE
    ind_break <- tryCatch(
      expr = {
        out=glmnet::glmnet(x=xmat,y=dlab,family="binomial",weights=weight,alpha=0,lambda=Lambda,standardize = FALSE)
      },
      error = function(e){
        print("ERROR, from SETMnull")
        return(TRUE)
      },
      warning = function(w){
        print(" Warning, from SETMnull")
        return(TRUE)
      })

    if(isTRUE(ind_break)){
      print("MY EXIT")
      return(list( penal_loglik=NA))
    }


    theta=as.numeric(coef(out))
    salp=theta[1]
    beta=theta[-1]


    hatpi1=1/(1+exp(salp+posdata%*%beta))
    hatpi2=1/(1+exp(salp+negdata%*%beta))
    hatpi=c(as.numeric(hatpi1),as.numeric(hatpi2))
    lam=sum(1-wi)/(n+sum(wi))
    n_hatpi=hatpi/sum(hatpi)
    alp=salp-log(lam)

    part1=sum(log(n_hatpi[1:n]))
    part2<-sum(as.numeric(log(pival/(1+exp(salp+negdata%*%beta))+
                                (1-pival)/(exp(-alp-negdata%*%beta)+lam)   )   ))

    newlik<-part1+part2-m*log(sum(hatpi))
    newlik_penal=newlik -(1/2*Lambda)*(beta%*%beta)
    err=newlik_penal-oldlik_penal

  }

  if(err<0){
    return(
      list(
        pival=old_pival,
        alp=old_alp,
        beta=old_beta,
        loglik=old_newlik,itenum=old_itenum,
        penal_loglik=oldlik_penal
      )
    )

  }
  else{
    return(
      list(
        pival=pival,
        alp=alp,
        beta=beta,
        loglik=newlik,itenum=itenum,
        penal_loglik=newlik_penal
      )
    )}

}
