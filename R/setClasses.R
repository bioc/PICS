#' @exportClass pics
setClass("pics", representation(estimates="list",score="numeric",scoreF="numeric",scoreR="numeric",Nmerged="numeric",converge="logical",range="numeric",chr="character"),
		prototype(estimates=list(w=numeric(0),mu=numeric(0),delta=numeric(0),sigmaSqF=numeric(0),sigmaSqR=numeric(0),seMu=numeric(0),seMuF=numeric(0),seMuR=numeric(0)),score=numeric(0),scoreF=numeric(0),scoreR=numeric(0),Nmerged=numeric(0),converge=logical(0),range=numeric(0),chr=character(0)))

#' @exportClass picsError
setClass("picsError", representation(errorCode="character"),prototype(errorCode=character(0)))

#' @exportClass picsList
setClass("picsList", representation(List="list", paraPrior="list", paraEM="list", minReads="list", N="integer", Nc="integer"), 
		prototype(List=list(0),minReads=list(perPeak=integer(0),perRegion=integer(0)), paraPrior=list(xi=double(0),rho=double(0),alpha=double(0),beta=double(0)),paraEM=list(kMax=integer(0),B=integer(0),tol=double(0)),N=integer(0), Nc=integer(0)))

### Constructor
#' @export
newPics<-function(w,mu,delta,sigmaSqF,sigmaSqR,seMu,seMuF,seMuR,score,scoreF,scoreR,Nmerged,converge,range,chr){
  if(!all(is.double(w))){
    stop("Argument 'w' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(mu))){
    stop("Argument 'mu' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(delta))){
    stop("Argument 'delta' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(sigmaSqF)) | !all(is.double(sigmaSqR))){
    stop("Argument 'sigmaSqF/sigmaSqR' must be numeric ", call.=FALSE)
  }
  if(!all(is.double(score))){
    stop("Argument 'score' must be numeric ", call.=FALSE)
  }
  if(!is.numeric(Nmerged)){
    stop("Argument 'Nmerged' must be numeric ", call.=FALSE)
  }
  if(!is.logical(converge)){
    stop("Argument 'converge' must be logical ", call.=FALSE)
  }
  if(!is.character(chr)){
    stop("Argument 'chr' must be a character string", call.=FALSE)
  }
  new("pics", estimates=list(w=w,mu=mu,delta=delta,sigmaSqF=sigmaSqF,sigmaSqR=sigmaSqR,seMu=seMu,seMuF=seMuF,seMuR=seMuR),converge=converge,score=score,scoreF=scoreF,scoreR=scoreR,Nmerged=Nmerged,range=range,chr=chr)
}

# In case the algorithm does not converge
#' @export
newPicsError<-function(string){
  if(!is.character(string)){
    stop("Argument 'errorCode' must be of class character", call.=FALSE)
  }
  new("picsError", errorCode=string)
}

#' @export
newPicsList<-function(List, paraEM, paraPrior, minReads, N, Nc){
  if(!is.list(paraEM) & !all(sapply(paraEM,"is.numeric"))){
    stop("Argument 'paraEM' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(paraPrior) & !all(sapply(paraPrior,"is.numeric"))){
    stop("Argument 'paraPrior' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(minReads) & !all(sapply(minReads,"is.numeric"))){
    stop("Argument 'minReads' must be a list of numeric arguments", call.=FALSE)
  }
  if(!all((lapply(List,"class")=="pics" | lapply(List,"class")=="picsError"))){
    stop("Argument 'List' must be a list of 'pics' or 'picsError' arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc)){
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  new("picsList", List=List, paraEM=paraEM, paraPrior=paraPrior, minReads=minReads, N=N, Nc=Nc)
}