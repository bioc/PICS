#' Classe and methods for candidate regions from paired-end data
#' 
#' A segReadsPE object represents a single candoidate region, including all its
#' informative reads and mappability profile.
#' 
#' @exportClass segReadsPE
setClass("segReadsPE", contains = "segReads", representation(yFm = "numeric", yRm = "numeric", cFm = "numeric", cRm = "numeric"), 
    prototype(yFm = numeric(0), yRm = numeric(0), cFm = numeric(0), cRm = numeric(0)))

#' Class and methods for list of candidate regions from paired-end data
#' 
#' @exportClass segReadsListPE
setClass("segReadsListPE", contains = "segReadsList", representation(NFm = "integer", NRm = "integer", NcFm = "integer", NcRm = "integer"), 
    prototype(list(List = list(0), paraSW = list(islandDepth = integer(0), min_cut = integer(0), max_cut = integer(0), xi = 0), 
        NFm = integer(0), NRm = integer(0), NcFm = integer(0), NcRm = integer(0))))

#' @describeIn segReadsPE-class \code{segReadsPE} constructor.
#' @export
segReadsPE <- function(yF, yR, yFm, yRm, cF, cR, cFm, cRm, map, chr) {
    if (!is.vector(yF) || !is.vector(yR) || !is.numeric(yF) || !is.numeric(yR)) {
        stop("Argument 'yF/yR' must be numeric vectors ", call. = FALSE)
    }
    if (!is.vector(yFm) || !is.vector(yRm) || !is.numeric(yFm) || !is.numeric(yRm)) {
        stop("Argument 'yFm/yRm' must be numeric vectors ", call. = FALSE)
    }
    if ((!is.vector(cF) || !is.vector(cR) || !is.numeric(cF) || !is.numeric(cR)) & (!is.null(cF) || !is.null(cR))) {
        stop("Argument 'cF/cR' must be numeric vectors ", call. = FALSE)
    }
    if ((!is.vector(cFm) || !is.vector(cRm) || !is.numeric(cFm) || !is.numeric(cRm)) & (!is.null(cFm) || !is.null(cRm))) {
        stop("Argument 'cFm/cRm' must be numeric vectors ", call. = FALSE)
    }
    if (!is.matrix(map)) {
        stop("Argument 'map' must be a matrix ", call. = FALSE)
    }
    new("segReadsPE", yF = yF, yR = yR, yFm = yFm, yRm = yRm, cF = cF, cR = cR, cFm = cFm, cRm = cRm, map = map, chr = chr)
}

#' @describeIn segReadsListPE-class \code{segReadsListPE} constructor.
#' @export
segReadsListPE <- function(List, paraSW, N, NFm, NRm, Nc, NcFm, NcRm) {
    if (!is.list(paraSW) & !all(sapply(paraSW, "is.numeric"))) {
        stop("Argument 'paraSW' must be a list of numeric arguments", call. = FALSE)
    }
    if (any(lapply(List, "class") != "segReadsPE")) {
        stop("Argument 'List' must be a list of segReadsPE arguments", call. = FALSE)
    }
    if (!is.integer(N) | !is.integer(Nc)) {
        stop("Argument 'N' and 'Nc' must be integers", call. = FALSE)
    }
    if (!is.integer(NFm) | !is.integer(NRm)) {
        stop("Argument 'NFm' and 'NRm' must be integers", call. = FALSE)
    }
    if (!is.integer(NcFm) | !is.integer(NcRm)) {
        stop("Argument 'NcFm' and 'NcRm' must be integers", call. = FALSE)
    }
    
    new("segReadsListPE", List = List, paraSW = paraSW, N = N, NFm = NFm, NRm = NRm, Nc = Nc, NcFm = NcFm, NcRm = NcRm)
}

## METHODS ##

#' @describeIn segReadsListPE-class subset method
#' @param x,i,j,...,drop,exact Arguments passed to subset methods
#' @export 
setMethod("[", "segReadsListPE", function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) {
        return(x)
    }
    if (!missing(j)) {
        stop("incorrect number of dimensions")
    } else {
        segReadsListPE(x@List[i], x@paraSW, x@N, x@NFm, x@NRm, x@Nc, x@NcFm, x@NcRm)
    }
})

#' @describeIn segReadsListPE-class subset method
#' @export 
setMethod("[[", "segReadsListPE", function(x, i, j, ..., exact = TRUE) {
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
