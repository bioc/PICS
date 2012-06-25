segmentPICS<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, dataType="TF",maxLregion=0,minLregion=100)
{
  #Paras depends on the datatype
  step=20
  if(dataType=="TF")  width=250
  if(dataType=="H")   width=150
  
  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=ReadsInRegion, jitter=FALSE, maxLregion=maxLregion,minLregion=minLregion,
			step=step, width=width, package="PICS")
  return(newSet)
}
