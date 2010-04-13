segmentReads<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3, jitter=FALSE, dataType="TF",maxLregion=0,minLregion=100)
{
  #maxLregion save max allowable region length, if it is not positive, that means no upper bound the regions length

  ## Check that we have the right data type
  if(!is(data,"GenomeData") | !is(data,"GenomeData"))
  {
    stop("The input data should be 'GenomeData' objects")
  }else
  {
  	data<-.switchReads(data)
  }
  ## Check that we have the same number of chromosomes
  if(length(data)!=length(dataC) & !is.null(dataC))
  {
    stop("Your IP and control data do not have the same number of chromosomes!")
  }
  if(is.null(dataC))
  {
    dataC<-GenomeData(vector("list",length(data)))
  }else
  {
  	dataC<-.switchReads(dataC)
  }
  
  if(dataType=="TF")  mywidth=250
  if(dataType=="H")   mywidth=150

  ## Total number of reads per sample
  lCont<-lapply(dataC,function(x){length(x[[1]])+length(x[[2]])})
  lIP<-lapply(data,function(x){length(x[[1]])+length(x[[2]])})

  if(is.null(minReads))
  {
    length<-lapply(data,function(x){diff(c(min(x[[1]][1],x[[2]][1]),max(tail(x[[1]],1),tail(x[[2]],1))))})
    minReads<-(ceiling(as.numeric(lIP)/as.numeric(length)*mywidth+1))
    minReads<-as.integer(names(which.max(table(minReads))))
    # print(minReads)
  }

  paraSW<-list(step=as.integer(20), width=as.integer(mywidth), minReads=as.integer(minReads))
  if(!is.null(map) & !is(map,"RangedData"))
  {
    stop("Map should be a 'RangedData' object")
  }
  else if(is.null(map))
  {
    start<-NULL
    end<-NULL
  }
  else
  {
  	tempS <- lapply(map,start)
  	tempE <- lapply(map,end)
  	start <- end <- vector("list",length(data))
  	chrs  <- names(start) <- names(end) <- names(data)
	for(cc in chrs)
	{
		start[[cc]]	<- tempS[[cc]]
		end[[cc]]	<- tempE[[cc]]
	}
  }
  
  #browser()
  if (maxLregion>0) maxStep=(maxLregion-2*paraSW$width)/paraSW$step else maxStep=0
  ## Perform the segmentation
  newSegReadsList<-.Call("segReadsAll", data, dataC, start, end, as.integer(jitter), paraSW , as.integer(maxStep), as.integer(minLregion),PACKAGE="PICS")


  temp<-unlist(newSegReadsList,recursive=FALSE,use.names=FALSE)
  if(is.null(temp))
  {
    stop("No Candidate regions found, you should decrease 'minReads'")
  }
  newSet<-segReadsList(temp,paraSW,as.integer(sum(unlist(lIP))),as.integer(sum(unlist(lCont))))
  
  ttt=.summarySeg(newSet)
  indrm=((ttt$L<minLregion)|(ttt$NF<minReadsInRegion)|(ttt$NR<minReadsInRegion))
  newSet@List=newSet@List[!indrm]
  
  return(newSet)
}

##check if the listData slot of the GenomeData 'x' has "-" before "+", if yes then switch
.switchReads <- function(x)
{
	if (names(x@listData[[1]])[[1]]=="-")
	{
		for(i in 1:length(x@listData))	x@listData[[i]]=x@listData[[i]][2:1]
	}
	return(x)
}

## summary a segmentList object return some information of each segment as a dataframe
## returned info including 
# chr: chromosome id
# NF : number of forward reads
# NR : number of reverse reads
# L  : length of segment
# min: start location of segments
# max: end location of segments
.summarySeg <- function(x)
{
	temp<-.Call("getSegL", x@List, PACKAGE="PICS");
	ans <- data.frame(chr=temp[[1]],NF=temp[[2]],NR=temp[[3]],L=temp[[4]],min=temp[[5]],max=temp[[6]])
	ans$chr <- as.character(ans$chr)
	return(ans)
}

