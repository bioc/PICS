#' Generics associated with pics classes
#' 
#' These generics are used in methods of \code{pics}, \code{picsError} and
#' \code{picsList} classes. See class help pages for method documentation.
#' 
#' @param x,... Object and argumemts passed to the methods.
#' 
#' @name pics-generics
NULL

#' @describeIn pics-generics Start accessor
#' @export
setGeneric("minRange", function(x, ...) standardGeneric("minRange"))

#' @describeIn pics-generics End accessor
#' @export
setGeneric("maxRange", function(x, ...) standardGeneric("maxRange"))

#' @describeIn pics-generics Score accessor
#' @export
setGeneric("score", function(x, ...) standardGeneric("score"))

#' @describeIn pics-generics Reverse score accessor
#' @export
setGeneric("scoreReverse", function(x, ...) standardGeneric("scoreReverse"))

#' @describeIn pics-generics Forward score accesor
#' @export
setGeneric("scoreForward", function(x, ...) standardGeneric("scoreForward"))

#' @describeIn pics-generics Chromosome accessor
#' @export
setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))

#' @describeIn pics-generics se accessor
#' @export
setGeneric("se", function(x, ...) standardGeneric("se"))

#' @describeIn pics-generics Forward se accessor
#' @export
setGeneric("seF", function(x, ...) standardGeneric("seF"))

#' @describeIn pics-generics Reverse se accessor
#' @export
setGeneric("seR", function(x, ...) standardGeneric("seR"))

#' @describeIn pics-generics sigmaSqF accessor
#' @export
setGeneric("sigmaSqF", function(x, ...) standardGeneric("sigmaSqF"))

#' @describeIn pics-generics sigmaSqR accessor
#' @export
setGeneric("sigmaSqR", function(x, ...) standardGeneric("sigmaSqR"))

#' @describeIn pics-generics delta accessor
#' @export
setGeneric("delta", function(x, ...) standardGeneric("delta"))

#' @describeIn pics-generics mu accessor
#' @export
setGeneric("mu", function(x, ...) standardGeneric("mu"))

#' @describeIn pics-generics w accessor
#' @export
setGeneric("w", function(x, ...) standardGeneric("w"))

#' @describeIn pics-generics K accessor
#' @export
setGeneric("K", function(x, ...) standardGeneric("K"))

#' @describeIn pics-generics Return error codes
#' @export
setGeneric("code", function(x, ...) standardGeneric("code"))


#' PICS class
#' 
#' This object is used to gather all parameters from fitting PICS to a single candidate region. 
#' 
#' @slot estimates A list containing all parameters estimates as well as standard errors.
#' @slot Nmerged The number of binding events that were merged; binding events that overlap are merged.
#' @slot converge A \code{logical} value indicating whether the EM algorithm has converged.
#' @slot chr The candidate region's chromosome.
#' @slot score Score of the binding event
#' @slot scoreF Forward score of the binding event.
#' @slot scoreR Reverse score of the binding event.
#' @slot range Genomic ranges.
#' 
#' @exportClass pics
setClass("pics", 
         representation(estimates = "list", score = "numeric", scoreF = "numeric", scoreR = "numeric", 
                        Nmerged = "numeric", converge = "logical", range = "numeric", chr = "character"),
         prototype(estimates = list(w = numeric(0), mu = numeric(0), delta = numeric(0),
                                    sigmaSqF = numeric(0), sigmaSqR = numeric(0), seMu = numeric(0), 
                                    seMuF = numeric(0), seMuR = numeric(0)),
                   score = numeric(0), scoreF = numeric(0), scoreR = numeric(0), 
                   Nmerged = numeric(0), converge = logical(0), range = numeric(0), 
                   chr = character(0)))

#' picsError class
#' 
#' This class is used in cases when the algorithm does not converge.
#' 
#' @slot errorCode The error code for debugging.
#' 
#' @exportClass picsError
setClass("picsError", representation(errorCode = "character"), prototype(errorCode = character(0)))

#' List of PICS objects
#' 
#' @exportClass picsList
setClass("picsList", representation(List = "list", paraPrior = "list", paraEM = "list", minReads = "list", N = "integer", Nc = "integer"), prototype(List = list(0), 
    minReads = list(perPeak = integer(0), perRegion = integer(0)), paraPrior = list(xi = double(0), rho = double(0), alpha = double(0), beta = double(0)), paraEM = list(kMax = integer(0), 
        B = integer(0), tol = double(0)), N = integer(0), Nc = integer(0)))

newPics <- function(w, mu, delta, sigmaSqF, sigmaSqR, seMu, seMuF, seMuR, score, scoreF, scoreR, Nmerged, converge, range, chr) {
    if (!all(is.double(w))) {
        stop("Argument 'w' must be numeric ", call. = FALSE)
    }
    if (!all(is.double(mu))) {
        stop("Argument 'mu' must be numeric ", call. = FALSE)
    }
    if (!all(is.double(delta))) {
        stop("Argument 'delta' must be numeric ", call. = FALSE)
    }
    if (!all(is.double(sigmaSqF)) | !all(is.double(sigmaSqR))) {
        stop("Argument 'sigmaSqF/sigmaSqR' must be numeric ", call. = FALSE)
    }
    if (!all(is.double(score))) {
        stop("Argument 'score' must be numeric ", call. = FALSE)
    }
    if (!is.numeric(Nmerged)) {
        stop("Argument 'Nmerged' must be numeric ", call. = FALSE)
    }
    if (!is.logical(converge)) {
        stop("Argument 'converge' must be logical ", call. = FALSE)
    }
    if (!is.character(chr)) {
        stop("Argument 'chr' must be a character string", call. = FALSE)
    }
    new("pics", estimates = list(w = w, mu = mu, delta = delta, sigmaSqF = sigmaSqF, sigmaSqR = sigmaSqR, seMu = seMu, seMuF = seMuF, seMuR = seMuR), converge = converge, 
        score = score, scoreF = scoreF, scoreR = scoreR, Nmerged = Nmerged, range = range, chr = chr)
}

newPicsError <- function(string) {
    if (!is.character(string)) {
        stop("Argument 'errorCode' must be of class character", call. = FALSE)
    }
    new("picsError", errorCode = string)
}

newPicsList <- function(List, paraEM, paraPrior, minReads, N, Nc) {
    if (!is.list(paraEM) & !all(sapply(paraEM, "is.numeric"))) {
        stop("Argument 'paraEM' must be a list of numeric arguments", call. = FALSE)
    }
    if (!is.list(paraPrior) & !all(sapply(paraPrior, "is.numeric"))) {
        stop("Argument 'paraPrior' must be a list of numeric arguments", call. = FALSE)
    }
    if (!is.list(minReads) & !all(sapply(minReads, "is.numeric"))) {
        stop("Argument 'minReads' must be a list of numeric arguments", call. = FALSE)
    }
    if (!all((lapply(List, "class") == "pics" | lapply(List, "class") == "picsError"))) {
        stop("Argument 'List' must be a list of 'pics' or 'picsError' arguments", call. = FALSE)
    }
    if (!is.integer(N) | !is.integer(Nc)) {
        stop("Argument 'N' and 'Nc' must be integers", call. = FALSE)
    }
    new("picsList", List = List, paraEM = paraEM, paraPrior = paraPrior, minReads = minReads, N = N, Nc = Nc)
}

## METHODS ##

#' @describeIn picsError-class show method
#' @param object,x A \code{picsError} object.
#' @export
setMethod("show", "picsError",
          function(object){
            cat("Object of class ",class(object),"\n")
            cat("This object has the following slot: \n")
            cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
          })

#' @describeIn pics-class show method
#' @param object,x A \code{pics} object.
#' @export
setMethod("show", "pics",
          function(object){
            cat("Object of class ",class(object),"\n")
            cat("This object has the following slots: \n")
            cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
            })

#' @describeIn picsList-class show method
#' @param object,x A \code{pics} object.
#' @export
setMethod("show", "picsList", function(object){
  cat("Object of class ",class(object),"\n")
  cat("This object has the following slots: \n")
  cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
  cat("List is a list of 'pics' or picsError ojects\n")
  })

#' @describeIn pics-class Get start of range
setMethod("minRange", "pics", function(x){return(x@range[1])})
#' @describeIn picsError-class Get start of range
setMethod("minRange", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Get start of range
setMethod("minRange", "picsList", function(x){
  ans <- .Call("getMin", x@List, PACKAGE="PICS");
  return(ans)
})

#' @describeIn pics-class Get end of range
setMethod("maxRange", "pics", function(x){return(x@range[2])})
#' @describeIn picsError-class Get end of range
setMethod("maxRange", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Get end of range
setMethod("maxRange", "picsList",function(x){
  ans<-.Call("getMax", x@List, PACKAGE="PICS");
  return(ans)
})

#' @describeIn pics-class Score accessor.
setMethod("score", "pics", function(x){return(x@score)})
#' @describeIn picsError-class Score accessor.
setMethod("score", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Score accessor.
setMethod("score", "picsList", function(x){
  ans<-.Call("getScore", x@List, PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics Score accessor
setMethod("score", "data.frame",function(x){return(x$score)})


#' @describeIn pics-class Reverse score accessor.
setMethod("scoreReverse", "pics",function(x){return(x@scoreR)})
#' @describeIn picsError-class Reverse score accessor.
setMethod("scoreReverse", "picsError",function(x){return(NULL)})
#' @describeIn picsList-class Reverse score accessor.
setMethod("scoreReverse", "picsList", function(x){
  ans<-.Call("getScoreR", x@List, PACKAGE="PICS");
  return(ans)            
})
#' @describeIn pics-generics Reverse score accessor
setMethod("scoreReverse", "data.frame", function(x){return(x$scoreR)})

#' @describeIn pics-class Forward score accessor.
setMethod("scoreForward", "pics", function(x){return(x@scoreF)})
#' @describeIn picsError-class Forward score accessor.
setMethod("scoreForward", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Forward score accessor.
setMethod("scoreForward", "picsList", function(x){
  ans<-.Call("getScoreF", x@List, PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics Forward score accessorse
setMethod("scoreForward", "data.frame", function(x){return(x$scoreF)})

#' @describeIn pics-class Chromosome accessor
setMethod("chromosome", "pics", function(x){return(rep(x@chr, length(x@estimates$w)))})
#' @describeIn picsError-class Chromosome accessor
setMethod("chromosome", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Chromosome accessor
setMethod("chromosome", "picsList", function(x){
  ans<-.Call("getChr", x@List, PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics chromosome accessor
setMethod("chromosome", "data.frame", function(x){return(x$chr)})

#' @describeIn pics-class se accessor
setMethod("se", "pics",function(x){return(x@estimates$seMu)})
#' @describeIn picsError-class se accessor
setMethod("se", "picsError",function(x){return(NULL)})
#' @describeIn picsList-class se accessor
setMethod("se", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(5), PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics se accessor
setMethod("se", "data.frame", function(x){return(x$se)})

#' @describeIn pics-class Forward se accessor
setMethod("seF", "pics", function(x){return(x@estimates$seMuF)})
#' @describeIn picsError-class Forward se accessor
setMethod("seF", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Forward se accessor
setMethod("seF", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(6), PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics seF accessor
setMethod("seF", "data.frame", function(x){return(x$seF)})

#' @describeIn pics-class Reverse se accessor
setMethod("seR", "pics", function(x){return(x@estimates$seMuR)})
#' @describeIn picsError-class Reverse se accessor
setMethod("seR", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class Reverse se accessor
setMethod("seR", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(7), PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics seR accessor
setMethod("seR", "data.frame", function(x){return(x$seR)})

#' @describeIn pics-class sigmaSqF accessor
setMethod("sigmaSqF", "pics", function(x){return(x@estimates$sigmaSqF)})
#' @describeIn picsError-class sigmaSqF accessor
setMethod("sigmaSqF", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class sigmaSqF accessor
setMethod("sigmaSqF", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(3), PACKAGE="PICS");
  return(ans)
})
#' @describeIn pics-generics sigmaSqF accessor
setMethod("sigmaSqF", "data.frame", function(x){return(x$sigmaSqF)})

#' @describeIn pics-class sigmaSqR accessor
setMethod("sigmaSqR", "pics", function(x){return(x@estimates$sigmaSqR)})
#' @describeIn picsError-class sigmaSqR accessor
setMethod("sigmaSqR", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class sigmaSqR accessor
setMethod("sigmaSqR", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(4), PACKAGE="PICS")
  return(ans)
})
#' @describeIn pics-generics sigmaSqR accessor
setMethod("sigmaSqR", "data.frame",	function(x){return(x$sigmaSqR)})

#' @describeIn pics-class delta accessor
setMethod("delta", "pics", function(x){return(x@estimates$delta)})
#' @describeIn picsError-class delta accessor
setMethod("delta", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class delta accessor
setMethod("delta", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(2), PACKAGE="PICS")
  return(ans)
})
#' @describeIn pics-generics delta accessor
setMethod("delta", "data.frame", function(x){return(x$delta)})


#' @describeIn pics-class mu accessor
setMethod("mu", "pics", function(x){return(x@estimates$mu)})
#' @describeIn picsError-class mu accessor
setMethod("mu", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class mu accessor
setMethod("mu", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(1), PACKAGE="PICS")
  return(ans)
})
#' @describeIn pics-generics mu accessor
setMethod("mu", "data.frame", function(x){return(x$mu)})

#' @describeIn pics-class w accessor
setMethod("w", "pics", function(x){return(x@estimates$w)})
#' @describeIn picsError-class w accessor
setMethod("w", "picsError", function(x){return(NULL)})
#' @describeIn picsList-class w accessor
setMethod("w", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(0), PACKAGE="PICS")
  return(ans)
})

#' @describeIn pics-class K accessor
setMethod("K", "pics", function(x){return(length(x@estimates$w))})
#' @describeIn picsError-class K accessor
setMethod("K", "picsError", function(x){return(0)})
#' @describeIn picsList-class K accessor
setMethod("K", "picsList", function(x){
  ans<-.Call("getK", x@List, PACKAGE="PICS")
  return(ans)
})

#' @describeIn pics-class Error code accessor
setMethod("code", "pics", function(x){return("")})
#' @describeIn picsError-class Error code accessor
setMethod("code", "picsError", function(x){return(x@errorCode)})
#' @describeIn picsList-class Error code accessor
setMethod("code", "picsList", function(x){
  temp<-lapply(x@List,"code")
  return(unlist(temp))
})

#' @describeIn picsList-class Return the length of the object
setMethod("length", "picsList", function(x){return(length(x@List))})

#' @describeIn picsList-class Summary of the object
setMethod("summary", "picsList", function(object) {
  cat("** Experiment information ** \n")
  cat("Chromosomes interogated: ")
  cat(unique(chromosome(object)), "\n")
  cat("Number of reads:")
  cat("In IP: ", object@N, " in control: ", object@Nc, "\n")
  cat("** Prior parameters ** \n")
  cat("The following settings were used:\n")
  cat("  Hyper parameters for the fragment length distribution:\n")
  cat("  xi, rho, alpha, beta: ", object@paraPrior$xi, ",", object@paraPrior$rho, ",", object@paraPrior$alpha, ",", object@paraPrior$beta, "\n")
  cat("** Score summary ** \n")
  print(summary(score(object)))
  cat("** Fragment length distribution summary ** \n")
  print(summary(delta(object)))
  cat("** Summary on the number of binding events per candidate region** \n")
  summary(K(object))
})

#' @describeIn pics-class Summary of the object
setMethod("summary", "pics", function(object) {
  cat("** Score ** \n")
  cat(score(object), "\n")
  cat("** Fragment length estimate ** \n")
  cat(delta(object), "\n")
  cat("** Number of binding events in the candidate region** \n")
  cat(K(object), "\n")
})

#' @describeIn picsList-class Subset list
#' @param i,j,...,drop,exact Arguments passed to subset functions
setMethod("[", "picsList", function(x, i, j, ..., drop = FALSE) {
  if (missing(i)) {
    return(x)
  }
  if (!missing(j)) {
    stop("incorrect number of dimensions")
  } else {
    newPicsList(x@List[i], x@paraEM, x@paraPrior, x@minReads, x@N, x@Nc)
  }
})

#' @describeIn picsList-class Subset element
setMethod("[[", "picsList", function(x, i, j, ..., exact = TRUE) {
  if (length(i) != 1) {
    stop("subscript out of bounds (index must have length 1)")
  }
  if (missing(i)) {
    return(x)
  }
  if (!missing(j)) {
    stop("incorrect number of dimensions")
  }
  x@List[[i]]
})