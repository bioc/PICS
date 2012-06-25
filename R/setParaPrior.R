# Default prior values for PICS
paraPriorTF<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=0)


setParaPrior<-function(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
{
  if(dataType!="TF")
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
  
  unlockBinding("paraPriorTF", environment(PICS))
  assign("paraPriorTF",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=environment(PICS))
  assign("paraPriorTF",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=.GlobalEnv)
  lockBinding("paraPriorTF", environment(PICS))
}
