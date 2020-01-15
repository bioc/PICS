#' Classes and functions to segment the genome in candidate regions
#' 
#' Pre-process bidirectional aligned reads data from a single ChIP-Seq 
#' experiment to detect candidate regions with a minimum number of forward and 
#' reverse reads. These candidate regions will then be processed by \code{PICS}.
#' 
#' @exportClass segReads
setClass("segReads", representation(yF = "numeric", yR = "numeric", cF = "numeric", cR = "numeric", map = "matrix", 
    chr = "character"), prototype(yF = numeric(0), yR = numeric(0), cF = numeric(0), cR = numeric(0), map = matrix(0, 
    0, 2), chr = character(0)))

#' @rdname segReads-class
#' @exportClass segReadsList
setClass("segReadsList", representation(List = "list", paraSW = "list", N = "integer", Nc = "integer"), prototype(List = list(List = list(0), 
    paraSW = list(step = integer(0), width = integer(0), minReads = integer(0)), N = integer(0), Nc = integer(0))))

#' @param yF A \code{numeric} vector. Forward reads.
#' @param yR A \code{numeric} vector. Reverse reads.
#' @param cF A \code{numeric} vector. Forward reads for the controls.
#' @param cR A \code{numeric} vector. Reverse reads for the controls.
#' @param map A \code{matrix}. The mappability profile.
#' @param chr A \code{character}. Chromosome name.
#' 
#' @describeIn segReads-class \code{segReads} Constructor
#' @export
segReads <- function(yF, yR, cF, cR, map, chr) {
    if (!is.vector(yF) || !is.vector(yR) || !is.numeric(yF) || !is.numeric(yR)) {
        stop("Argument 'yF/yR' must be numeric vectors ", call. = FALSE)
    }
    if ((!is.vector(cF) || !is.vector(cR) || !is.numeric(cF) || !is.numeric(cR)) & (!is.null(cF) || !is.null(cR))) {
        stop("Argument 'cF/cR' must be numeric vectors ", call. = FALSE)
    }
    
    if (!is.matrix(map)) {
        stop("Argument 'map' must be a matrix ", call. = FALSE)
    }
    if (!is.character(chr)) {
        stop("Argument 'chr' must be a character string", call. = FALSE)
    }
    
    new("segReads", yF = yF, yR = yR, cF = cF, cR = cR, map = map, chr = chr)
}

#' @param List A \code{list} of \code{segReads} objects.
#' @param paraSW A \code{list} of parameters for the genomic regions.
#' @param N A \code{numeric}. The number of reads in the data.
#' @param Nc A \code{numeric}. The number of reads in the control.
#' 
#' @note segReads and segReadsList objects are not meant to be built via the
#' constructors. The constructors are used in \code{segmentPICS}.
#' 
#' @describeIn segReads-class \code{segReadsList} Constructor
#' @export
segReadsList <- function(List, paraSW, N, Nc) {
    if (!is.list(paraSW) & !all(sapply(paraSW, "is.numeric"))) {
        stop("Argument 'paraSW' must be a list of numeric arguments", call. = FALSE)
    }
    if (any(lapply(List, "class") != "segReads")) {
        stop("Argument 'List' must be a list of segReads arguments", call. = FALSE)
    }
    if (!is.integer(N) | !is.integer(Nc)) {
        stop("Argument 'N' and 'Nc' must be integers", call. = FALSE)
    }
    new("segReadsList", List = List, paraSW = paraSW, N = N, Nc = Nc)
}


## METHODS ##

#' @describeIn segReads-class show method
#' @param object,x A \code{segReads} object.
#' @export
setMethod("show", "segReads", function(object) {
    cat("Object of class ", as.character(class(object)), "\n")
    cat("This object has the following slots: \n")
    cat(paste(names(getSlots(class(object))), collapse = ", "), "\n")
})

#' @describeIn segReads-class show method
#' @export
setMethod("show", "segReadsList", function(object) {
    cat("Object of class", as.character(class(object)), "\n")
    cat("This object has the following slots: \n")
    cat(paste(names(getSlots(class(object))), collapse = ", "), "\n")
    cat("List is a list of 'segReads' ojects, each of which has the following slots:\n")
    cat("yR, yF, cR, cF, map, chr\n")
})

#' @describeIn segReads-class map generic
setGeneric("map", function(x, ...) standardGeneric("map"))
#' @describeIn segReads-class map method
#' @export
setMethod("map", "segReads", function(x) {
    if (is.null(x@map) | (nrow(x@map) == 0)) {
        return(0)
    } else {
        n <- nrow(x@map)
        m <- min(x@yF[1], x@yR[1], x@map[1, 1])
        M <- max(tail(x@yF, 1), tail(x@yR, 1), x@map[n, 2])
        return(sum(diff(t(x@map)))/max(M - m, 1))
    }
})
#' @describeIn segReads-class map method
#' @export
setMethod("map", "segReadsList", function(x) {
    ans <- .Call("getMap", x@List, PACKAGE = "PICS")
    return(ans)
})

#' @describeIn segReads-class Return length of segReadsList
setMethod("length", "segReadsList", function(x){return(length(x@List))})

#' @describeIn segReads-class Summary method
#' @export
setMethod("summary", "segReadsList", function(object){
  cat("** Experiment information ** \n")
  cat("Chromosomes interogated: ")
  cat(unique(unlist(lapply(object@List,function(obj){obj@chr}))),"\n")
  cat("Number of reads")
  cat(" in IP: ",object@N," and in control: ",object@Nc,"\n")
  cat("** Segmentation parameters ** \n")
  cat("The following settings were used:\n")          
  cat("  Sliding window half width: ", object@paraSW$width,"\n")
  cat("  Step size: ", object@paraSW$step,"\n")          
  cat("  Minimum number of reads: ", object@paraSW$minReads,"\n")
  cat("** Segmentation summary ** \n")                    
  cat("Number of segmented regions:",length(object@List),"\n")
  cat("Summary on the number of Forward/Reverse reads per region:\n")
  cat("  Forward:\n") 
  cat("  ")
  tempF<-lapply(object@List,function(obj){length(obj@yF)})
  print(summary(as.integer(unlist(tempF))))
  cat("  Reverse:\n") 
  cat("  ")
  tempR<-lapply(object@List,function(obj){length(obj@yR)})
  print(summary(as.integer(unlist(tempR))))
  cat("Summary on the number of control Forward/Reverse reads per region:\n")
  cat("  Forward:\n") 
  cat("  ")
  tempF<-lapply(object@List,function(obj){length(obj@cF)})
  print(summary(as.integer(unlist(tempF))))
  cat("  Reverse:\n") 
  cat("  ")
  tempR<-lapply(object@List,function(obj){length(obj@cR)})
  print(summary(as.integer(unlist(tempR))))                    
  tempMap<-map(object)
  cat("** Mappability summary **\n")
  cat("Non mappable intervals cover an average ", mean(unlist(tempMap)),"% of all regions \n")                  
})

#' @describeIn segReads-class Summary method
#' @export
setMethod("summary", "segReads",
          function(object)
          {
            m<-min(object@yF[1],object@yR[1])
            M<-max(tail(object@yF,1),tail(object@yR,1))
            cat("** Region summary ** \n")
            cat("Summary on Forward reads:\n")
            print(summary(object@yF,digits=100))
            cat("Summary on Reverse reads:\n")
            print(summary(object@yR,digits=100))
            cat("Summary on control Forward reads:\n")
            print(summary(object@cF,digits=100))
            cat("Summary on control Reverse reads:\n")
            print(summary(object@cR,digits=100))
            cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M-m),"% of the region \n")
          })


#' @describeIn segReads-class Subset methods
#' @param i,j,...,exact,drop Additional arguments passed to subset methods.
#' @export
setMethod("[","segReadsList", function(x,i, j,..., drop=FALSE){
    if(missing(i)){
        return(x)
    }
    if(!missing(j)){
        stop("incorrect number of dimensions")
    } else {
        segReadsList(x@List[i],x@paraSW,x@N,x@Nc)
    }
})

#' @describeIn segReads-class Subset methods
#' @export
setMethod("[[","segReadsList", function(x, i, j, ..., exact = TRUE){
    if(length(i) != 1){
        stop("subscript out of bounds (index must have length 1)")
    }
    if(missing(i)){
        return(x)
    }
    if(!missing(j)){
        stop("incorrect number of dimensions")
    }
    x@List[[i]]
})


