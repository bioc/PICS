# Default prior values for PICS
paraEMTF<-list(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE)
paraEMH<-list(minK=1,maxK=15,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=FALSE,mapCorrect=TRUE)


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
  
  if(dataType=="TF")
  {
  unlockBinding("paraEMTF", environment(PICS))
  assign("paraEMTF",list(minK=minK,maxK=maxK,tol=tol,B=B,mSelect=mSelect,mergePeaks=mergePeaks,mapCorrect=mapCorrect),envir=environment(PICS))
  lockBinding("paraEMTF", environment(PICS))
  }
  else if(dataType=="H")
  {
    unlockBinding("paraEMH", environment(PICS))
    assign("paraEMH",list(minK=minK,maxK=maxK,tol=tol,B=B,mSelect=mSelect,mergePeaks=mergePeaks,mapCorrect=mapCorrect),envir=environment(PICS))
    lockBinding("paraEMH", environment(PICS))
  }
}
