## Functions to return a list of parameters to be used by PICS functions

setParaEM<-function(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType="TF")
{
  if(dataType!="TF" & dataType!="H")
  {
    stop("Object 'dataType' must be either 'TF' or 'H'")
  }
  if(!is.finite(tol) & tol<=0 & tol>1)
  {
    stop("'tol' must be a positive number between 0 and 1")
  }
  if(!is.finite(maxK) & maxK<=0 & round(maxK)!=maxK)
  {
    stop("'maxK' must be a positive integer")
  }
  if(!is.finite(B) & B<1 & round(B)!=B)
  {
    stop("'B' must be a positive integer")
  }
  if(!is.character(mSelect) & ((mSelect!="BIC") | (mSelect!="AIC3")))
  {
    stop("'lambda' must be a positive number")
  }
  if(!is.logical(mergePeaks))
  {
    stop("'mergePeaks' must be a logical value")
  }
  if(!is.logical(mapCorrect))
  {
    stop("'mapCorrect' must be a logical value")
  }
  return(list(minK=minK,maxK=maxK,tol=tol,B=B,mSelect=mSelect,mergePeaks=mergePeaks,mapCorrect=mapCorrect))
}

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
  return(list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu))
}

