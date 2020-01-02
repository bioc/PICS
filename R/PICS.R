#' @useDynLib PICS, .registration = TRUE
#' @import methods
NULL

#' Estimation of binding site positions
#' 
#' This object contains Estimation of binding site positions and has the
#' following slots: segReadsList, dataType.
#' 
#' @param segReadsList This object contains segmentation of Genome
#' @param dataType A \code{character}. The type of data you are processing: 
#'  'TF' for transcription factor.
#' @param paraEM A \code{list} of parameters for the EM algorithm as returned by 
#'  the \code{setParaEM} function. The default parameters should be good enough 
#'  for most usages.
#' @param paraPrior A \code{list} of parameters for the priors as returned by
#'  \code{setParaPrior}.
#' @param nCores An \code{integer}. The number of cores that should be used in
#'  parallel by the function.
#'
#' @return An object of class \code{picsList} containing the estimated binding
#' site positions.
#' 
#' @export
PICS <- function(segReadsList, dataType = NULL, paraEM = NULL, paraPrior = NULL, nCores = 1) {
    ### Constant used in the calculations
    cst <- gamma(3.5)/gamma(3)/sqrt(pi)
    minReads <- list(perPeak = 3, perRegion = 4)
    if (length(paraEM) != 7) {
        paraEM <- setParaEM(dataType = dataType)  #using PICS default paraEM
    }
    if (length(paraPrior) != 6) {
        paraPrior <- setParaPrior(dataType = dataType)  #using PICS default paraPrior
    }
    
    if (nCores > 1 & requireNamespace("parallel", quietly = TRUE)) {
        availCores <- parallel::detectCores()
        if (nCores > availCores) {
            warning("The number of cores required is higher than the available cores on this machine (", availCores, ").\n", immediate. = TRUE)
            nCores <- availCores
        }
        message("Using the parallel version of PICS with ", nCores, " cpus or cores")
        # Split into nCores segReadsList
        cl <- parallel::makeCluster(getOption("cl.cores", nCores))
        segSplit <- split(segReadsList, cut(1:length(segReadsList), nCores))
        # Use parallel version of lapply
        res <- unlist(parallel::parLapply(cl, segSplit, .fitModelAllkSplit, paraEM, paraPrior, minReads), recursive = FALSE)
        parallel::stopCluster(cl)
    } else {
        message("Using the serial version of PICS")
        res <- .Call("fitPICS", segReadsList, paraEM, paraPrior, minReads, PACKAGE = "PICS")
    }
    
    myPicsList <- newPicsList(res, paraEM, paraPrior, minReads, segReadsList@N, segReadsList@Nc)
    return(myPicsList)
}

.fitModelAllkSplit <- function(segReadsList, paraEM, paraPrior, minReads) {
    res <- .Call("fitPICS", segReadsList, paraEM, paraPrior, minReads, PACKAGE = "PICS")
}


## This function could be used to simulate random reads in the case there are no background reads
backgroundSim <- function(dataF, dataR, mapPro = NULL, gapPro = NULL, pRetain = 0.01) {
    obj <- .C("backgroundSim", dataF = as.double(dataF), dataR = as.double(dataR), nF = as.integer(length(dataF)), nR = as.integer(length(dataR)), 
        as.integer(mapPro[, 1]), as.integer(mapPro[, 2]), as.integer(nrow(mapPro)), as.integer(gapPro[, 1]), as.integer(gapPro[, 
            2]), as.integer(nrow(gapPro)), as.double(pRetain), PACKAGE = "PICS")
    
    list(dataF = obj$dataF, dataR = obj$dataR)
}

# it filter the data.frame converted from pics object
.filterPICS <- function(ss, filter = list(delta = c(50, 250), sigmaSq = 22500, se = 50, mu = c(0, Inf), chr = NULL)) {
    ind1 <- (ss$delta >= filter$delta[1]) & (ss$delta <= filter$delta[2])
    ind2 <- (ss$sigmaSqF < filter$sigmaSq) & (ss$sigmaSqR < filter$sigmaSq)
    ind3 <- (ss$mu >= filter$mu[1]) & (ss$mu < filter$mu[2])
    ind4 <- (ss$se < filter$se)
    ind4[is.na(ind4)] <- TRUE  # do not filter by SE if it is not calculable
    ind5 <- rep(TRUE, nrow(ss))
    if (length(filter$chr) > 0) 
        ind5 <- (ss$chr %in% filter$chr)
    ans <- ss[ind1 & ind2 & ind3 & ind4 & ind5, ]
    return(ans)
}

# filter nucleosomes predicted outside of segment range.
.filterPICS2 <- function(ss) {
    ind1 <- (ss$mu <= ss$maxRange) & (ss$mu >= ss$minRange)
    ans <- ss[ind1, ]
    return(ans)
}
