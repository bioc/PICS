#' Classe and methods for candidate regions
#' 
#' A segReads object represents a single candoidate region, including all its
#' informative reads and mappability profile.
#' 
#' @exportClass segReads
setClass("segReads", representation(yF = "numeric", yR = "numeric", cF = "numeric", cR = "numeric", map = "matrix", chr = "character"), 
    prototype(yF = numeric(0), yR = numeric(0), cF = numeric(0), cR = numeric(0), map = matrix(0, 0, 2), chr = character(0)))

#' Class and methods for list of candidate regions
#' 
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
#' @note segReads objects are not meant to be built via the
#' constructors. The constructors are used in \code{segmentPICS}.
#' 
#' @describeIn segReads-class \code{segReads} constructor
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
#' @note segReadsList objects are not meant to be built via the
#' constructors. The constructors are used in \code{segmentPICS}.
#' 
#' @describeIn segReadsList-class \code{segReadsList} constructor
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
#' @param object A \code{segReads} object.
#' @export
setMethod("show", "segReads", function(object) {
    cat("Object of class ", as.character(class(object)), "\n")
    cat("This object has the following slots: \n")
    cat(paste(names(getSlots(class(object))), collapse = ", "), "\n")
})

#' @describeIn segReadsList-class show method
#' @param object A \code{segReadsList} object.
#' @export
setMethod("show", "segReadsList", function(object) {
    cat("Object of class", as.character(class(object)), "\n")
    cat("This object has the following slots: \n")
    cat(paste(names(getSlots(class(object))), collapse = ", "), "\n")
    cat("List is a list of 'segReads' ojects, each of which has the following slots:\n")
    cat("yR, yF, cR, cF, map, chr\n")
})


setGeneric("map", function(object, ...) standardGeneric("map"))
#' @describeIn segReads-class map method
#' @export
setMethod("map", "segReads", function(object) {
    if (is.null(object@map) | (nrow(object@map) == 0)) {
        return(0)
    } else {
        n <- nrow(object@map)
        m <- min(object@yF[1], object@yR[1], object@map[1, 1])
        M <- max(tail(object@yF, 1), tail(object@yR, 1), object@map[n, 2])
        return(sum(diff(t(object@map)))/max(M - m, 1))
    }
})
#' @describeIn segReadsList-class map method
#' @export
setMethod("map", "segReadsList", function(object) {
    ans <- .Call("getMap", object@List, PACKAGE = "PICS")
    return(ans)
})

setMethod("length", "segReadsList", function(x) {
    return(length(x@List))
})

#' @describeIn segReadsList-class summary method
#' @export
setMethod("summary", "segReadsList", function(object) {
    cat("** Experiment information ** \n")
    cat("Chromosomes interogated: ")
    cat(unique(unlist(lapply(object@List, function(obj) {
        obj@chr
    }))), "\n")
    cat("Number of reads")
    cat(" in IP: ", object@N, " and in control: ", object@Nc, "\n")
    cat("** Segmentation parameters ** \n")
    cat("The following settings were used:\n")
    cat("  Sliding window half width: ", object@paraSW$width, "\n")
    cat("  Step size: ", object@paraSW$step, "\n")
    cat("  Minimum number of reads: ", object@paraSW$minReads, "\n")
    cat("** Segmentation summary ** \n")
    cat("Number of segmented regions:", length(object@List), "\n")
    cat("Summary on the number of Forward/Reverse reads per region:\n")
    cat("  Forward:\n")
    cat("  ")
    tempF <- lapply(object@List, function(obj) {
        length(obj@yF)
    })
    print(summary(as.integer(unlist(tempF))))
    cat("  Reverse:\n")
    cat("  ")
    tempR <- lapply(object@List, function(obj) {
        length(obj@yR)
    })
    print(summary(as.integer(unlist(tempR))))
    cat("Summary on the number of control Forward/Reverse reads per region:\n")
    cat("  Forward:\n")
    cat("  ")
    tempF <- lapply(object@List, function(obj) {
        length(obj@cF)
    })
    print(summary(as.integer(unlist(tempF))))
    cat("  Reverse:\n")
    cat("  ")
    tempR <- lapply(object@List, function(obj) {
        length(obj@cR)
    })
    print(summary(as.integer(unlist(tempR))))
    tempMap <- map(object)
    cat("** Mappability summary **\n")
    cat("Non mappable intervals cover an average ", mean(unlist(tempMap)), "% of all regions \n")
})

#' @describeIn segReads-class summary method
#' @export
setMethod("summary", "segReads", function(object) {
    m <- min(object@yF[1], object@yR[1])
    M <- max(tail(object@yF, 1), tail(object@yR, 1))
    cat("** Region summary ** \n")
    cat("Summary on Forward reads:\n")
    print(summary(object@yF, digits = 100))
    cat("Summary on Reverse reads:\n")
    print(summary(object@yR, digits = 100))
    cat("Summary on control Forward reads:\n")
    print(summary(object@cF, digits = 100))
    cat("Summary on control Reverse reads:\n")
    print(summary(object@cR, digits = 100))
    cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M - m), "% of the region \n")
})


#' @describeIn segReadsList-class Subset methods
#' @param x,i,j,...,drop,exact Arguments passed to subset methods
#' @export
setMethod("[", "segReadsList", function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) {
        return(x)
    }
    if (!missing(j)) {
        stop("incorrect number of dimensions")
    } else {
        segReadsList(x@List[i], x@paraSW, x@N, x@Nc)
    }
})

#' @describeIn segReadsList-class Subset methods
#' @export
setMethod("[[", "segReadsList", function(x, i, j, exact = TRUE) {
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


