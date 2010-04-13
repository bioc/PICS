# Default prior values for PICS
#priottype=1 means use the new prior regulate both forward/reverse peaks, otherwise old prior is used to regulate sum of precision of peaks
#paraPriorH<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=-0.00001,dMu=200,priortype=1)
#paraPriorTF<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=0,priortype=1)
#I changed the value hyperparameter, 
# If both forward/reverse precision become centered at previous center of sum precision, 
# then rho should change from 1 to 0.5 to keep variance of delta not change.
# The previous variance of precision was too small, I change alpha, beta to make it value of sigma_f and sigma_r center @ 50, and CI= c(1,100)
paraPriorH<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=-0.000064,dMu=200)
paraPriorTF<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=0)


setParaPrior<-function(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
{
  if(dataType!="TF" & dataType!="H")
  {
    stop("Object 'dataType' must be either 'TF' or 'H'")
  }
  if(!is.finite(xi))
  {
    stop("'xi' must be a numeric value")
  }
  if(!is.finite(rho) & rho<=0)
  {
    stop("'rho' must be a positive number")
  }
  if(!is.finite(alpha) & alpha<=0)
  {
    stop("'alpha' must be a positive number")
  }
  if(!is.finite(beta) & beta<=0)
  {
    stop("'beta' must be a positive number")
  }
  if(!is.finite(lambda) & lambda<=0)
  {
    stop("'lambda' must be a positive number")
  }
  if(!is.finite(dMu) & dMu<=0)
  {
    stop("'dMu' must be a positive number")
  }
  
  if(dataType=="TF")
  {
  unlockBinding("paraPriorTF", environment(PICS))
  assign("paraPriorTF",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=environment(PICS))
  lockBinding("paraPriorTF", environment(PICS))
  }
  else if(dataType=="H")
  {
    if(!(lambda<0))
    {
      warning("You have selected the Histone option but lambda=0")
    }
    unlockBinding("paraPriorH", environment(PICS))
    assign("paraPriorH",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=environment(PICS))
    lockBinding("paraPriorH", environment(PICS))
  }
}
