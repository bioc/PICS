PICS<-function(segReadsList,dataType="TF")
{
  ### Constant used in the calculations
  cst<-gamma(3.5)/gamma(3)/sqrt(pi)
  minReads<-list(perPeak=3,perRegion=4)

  if(dataType!="TF")
  {
    stop("Object 'dataType' must be either 'TF'")
  }
  else
  {
    paraPrior<-paraPriorTF
    paraEM<-paraEMTF
  }
  
  if(length(grep("snowfall",loadedNamespaces()))==0 || !sfParallel())
  {
    # C version
    res<-.Call("fitPICS", segReadsList, paraEM, paraPrior, minReads, PACKAGE="PICS")
  }
  else
  {
    # Number of clusters
    nClust<-sfCpus()
    # Split into nClust segReadsList
    segSplit<-split(segReadsList,cut(1:length(segReadsList),nClust))
    # Use a parallel version
    res<-unlist(sfLapply(segSplit,.fitModelAllkSplit,paraEM,paraPrior,minReads),recursive=FALSE)
  }

  myPicsList<-newPicsList(res,paraEM,paraPrior,minReads,segReadsList@N,segReadsList@Nc)
  return(myPicsList)
}

.fitModelAllkSplit<-function(segReadsList,paraEM,paraPrior,minReads)
{
  res<-.Call("fitPICS", segReadsList, paraEM, paraPrior, minReads, PACKAGE="PICS")
}


## This function could be used to simulate random reads in the case there are no background reads
.background<-function(dataF, dataR, mapPro=NULL,gapPro=NULL,pRetain=0.01)
{
  obj<-.C("background",
  dataF=as.double(dataF),
  dataR=as.double(dataR),
  nF=as.integer(length(dataF)),
  nR=as.integer(length(dataR)),
  as.integer(mapPro[,1]),
  as.integer(mapPro[,2]),
  as.integer(nrow(mapPro)),
  as.integer(gapPro[,1]),
  as.integer(gapPro[,2]),
  as.integer(nrow(gapPro)),
  as.double(pRetain),
  PACKAGE="PICS")
  
  list(dataF=obj$dataF,dataR=obj$dataR)
}

#it filter the data.frame converted from pics object
.filterPICS <- function(ss,filter=list(delta=c(50,250),sigmaSq=22500, se=50, mu=c(0,Inf), chr=NULL))
{
	ind1	<- (ss$delta>=filter$delta[1])&(ss$delta<=filter$delta[2])
	ind2	<- (ss$sigmaSqF<filter$sigmaSq)&(ss$sigmaSqR<filter$sigmaSq)
	ind3	<- (ss$mu>=filter$mu[1])&(ss$mu<filter$mu[2])
	ind4	<- (ss$se<filter$se)
	ind4[is.na(ind4)]	<- TRUE  # do not filter by SE if it is not calculatable
	ind5	<- rep(TRUE,nrow(ss))
	if (length(filter$chr)>0) ind5 <- (ss$chr %in% filter$chr)
	ans		<- ss[ind1&ind2&ind3&ind4&ind5,]
	return(ans)
}

#filter nucleosomes predicted outside of segment range.
.filterPICS2 <- function(ss)
{
	ind1	<- (ss$mu<=ss$maxRange)&(ss$mu>=ss$minRange)
	ans		<- ss[ind1,]
	return(ans)
}
