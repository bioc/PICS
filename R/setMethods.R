#' @importFrom graphics abline axis grid layout lines mtext par plot points segments stripchart symbols title 
#' @importFrom utils tail
#' @importFrom stats density dt reshape

setAs("picsList", "GRanges", function(from) {
    makeGRangesOutput(from, type = "bed", filter = list(delta = c(50, 300), se = c(0, 50), sigmaSqF = c(0, 22500), sigmaSqR = c(0, 22500), score = c(1, Inf)), length = 100)
})

# Transform a .bed like data.frame into a GRanges object
setAs("data.frame", "GRanges", function(from) {
    if (length(from) < 4) 
        stra <- "*" else stra <- from[, 4]
    GRanges(ranges = IRanges(start = from[, 2], end = from[, 3]), seqnames = from[, 1], strand = stra)
})


## REMOVED (Maybe temporarily)

# @describeIn pics-generics Return density
# @export
setGeneric("wigDensity", function(x, ...) standardGeneric("wigDensity"))
# @describeIn pics-class Density
setMethod("wigDensity", "pics", function(x, strand = "+", step = 10, sum = FALSE, filter = NULL, scale = TRUE) {
  # Check that all filters are passed
  missingNames <- !c("delta", "sigmaSqF", "sigmaSqR", "se", "seF", "seR", "score") %in% names(filter)
  filter[c("delta", "sigmaSqF", "sigmaSqR", "se", "seF", "seR", "score")[missingNames]] <- list(c(0, Inf))
  if (strand == "+") {
    strand <- 1
  } else if (strand == "-") {
    strand <- -1
  } else if (strand == "*") {
    strand <- 0
  } else {
    stop("Strand must be either '+', '-' or '*'")
  }
  strand <- as.double(paste(strand, "1", sep = ""))
  ans <- .Call("getDensity", x, strand, step, filter, sum, scale, PACKAGE = "PICS")
  return(ans)
})
# @describeIn picsList-class Density
setMethod("wigDensity", "picsList", function(x, strand = "+", step = 10, sum = FALSE, filter = NULL, scale = TRUE) {
  # Check that all filters are passed
  missingNames <- !c("delta", "sigmaSqF", "sigmaSqR", "se", "seF", "seR", "score") %in% names(filter)
  filter[c("delta", "sigmaSqF", "sigmaSqR", "se", "seF", "seR", "score")[missingNames]] <- list(c(0, Inf))
  
  if (strand == "+") {
    strand <- 1
  } else if (strand == "-") {
    strand <- -1
  } else if (strand == "*") {
    strand <- 0
  } else {
    stop("Strand must be either '+', '-' or '*'")
  }
  ans <- .Call("getDensityList", x, strand, step, filter, sum, scale, PACKAGE = "PICS")
  return(ans)
})
# @describeIn picsError-class Density
setMethod("wigDensity", "picsError", function(x, strand = NULL, step = NULL, sum = NULL, filter = NULL) {
  return(NULL)
})


# I need a concatenate method to add pics objects


setAs("picsList", "data.frame", function(from) {
    ans <- data.frame(ID = rep(1:length(from), K(from)), chr = chromosome(from), w = w(from), mu = mu(from), delta = delta(from), sigmaSqF = sigmaSqF(from), sigmaSqR = sigmaSqR(from), 
        se = se(from), score = score(from), scoreF = scoreForward(from), scoreR = scoreReverse(from), minRange = minRange(from), maxRange = maxRange(from))
    ans$chr <- as.character(ans$chr)
    ans <- ans[is.finite(ans$mu), ]
    return(ans)
})

#' Plot methods for PICS objects
#' 
#' Methods to plot \code{pics} and \code{segReads} objects and derived classes.
#' 
#' @param x,y Objects
#' @param addKernel Add kernel density estimate to the plot.
#' @param addNucleosome Add a nucleosome track to the plot.
#' @param addSe Add standard error to the plot.
#' @param main Main title.
#' @param regionIndex Add region index to the plot.

#' @param ... Arguments to be passed to the plot method.
#' 
#' @name plot-pics
NULL

#' @describeIn plot-pics Plot method for \code{pics} and \code{segReads}
#' @importFrom grDevices grey
#' @export
setMethod("plot", signature("pics", "segReads"), function(x, y, addKernel = FALSE, addNucleosome = FALSE, addSe = TRUE, main = NULL, ...) {
    # Set outer and figure margins to reduce gap between plots
    if (addNucleosome) {
        nG <- 4
    } else {
        nG <- 2
    }
    par(oma = c(2.5, 5, 5, 5), mar = c(0, 5, 0, 0), cex.lab = 2)
    layout(matrix(1:nG, ncol = 1), heights = c(0.5, 0.2, 0.1, 0.1, 0.1))
    
    step <- 5
    .densityMix <- function(x, para) {
        v <- 4
        w <- para$w
        mu <- para$mu
        sigmaSq <- para$sigmaSq
        sigma <- sqrt(sigmaSq)
        xNorm <- outer(-mu, x, "+")/sigma
        return(colSums(w * dt(xNorm, df = v)/sigma))
    }
    
    yF <- y@yF
    yR <- y@yR
    cF <- y@cF
    cR <- y@cR
    map <- y@map
    m <- min(yF[1], yR[1]) - 100
    M <- max(tail(yF, 1), tail(yR, 1)) + 100
    
    paraR <- list(w = x@estimates$w, mu = x@estimates$mu + x@estimates$delta/2, sigmaSq = x@estimates$sigmaSqR)
    paraF <- list(w = x@estimates$w, mu = x@estimates$mu - x@estimates$delta/2, sigmaSq = x@estimates$sigmaSqF)
    
    dR <- .densityMix(seq(m, M, step), paraR)
    dF <- .densityMix(seq(m, M, step), paraF)
    maxRange <- max(c(dF, dR))
    plot(seq(m, M, step), dF, xlim = c(m, M), ylim = c(0, maxRange), lty = 2, type = "l", xlab = "", ylab = "density", xaxt = "n", axes = FALSE)
    title(main = main, outer = TRUE, cex.main = 2)
    axis(2)
    axis(1)
    
    
    lines(seq(m, M, step), dR, lty = 2, col = 2)
    
    # if(length(map)>0) { nMap<-nrow(map) for(i in 1:nMap) { segments(map[i,1], 0, map[i,2], 0,lwd=3,col=3) } }
    
    # Add kernel density estimate
    if ((addKernel == TRUE) & (length(yF) > 1 & length(yR) > 1)) {
        dkF <- density(yF)
        dkR <- density(yR)
        lines(dkF, lty = 3)
        lines(dkR, col = 2, lty = 3)
    }
    
    # Add single components and se's
    K <- length(x@estimates$w)
    for (k in 1:K) {
        paraR <- list(w = x@estimates$w[k], mu = x@estimates$mu[k] + x@estimates$delta[k]/2, sigmaSq = x@estimates$sigmaSqR[k])
        paraF <- list(w = x@estimates$w[k], mu = x@estimates$mu[k] - x@estimates$delta[k]/2, sigmaSq = x@estimates$sigmaSqF[k])
        
        dsR <- .densityMix(seq(m, M, step), paraR)
        dsF <- .densityMix(seq(m, M, step), paraF)
        
        lines(seq(m, M, step), dsF, lty = 1)
        lines(seq(m, M, step), dsR, col = 2, lty = 1)
    }
    
    stripchart(yF[1], pch = ">", method = "overplot", cex = 2, at = 0.55, add = FALSE, axes = FALSE, xlim = c(m, M), ylim = c(0, 1))
    if (length(map) > 0) {
        nMap <- nrow(map)
        symbols((map[, 1] + map[, 2])/2, rep(0.35, nMap), rectangles = cbind(map[, 2] - map[, 1], rep(0.6, nMap)), inches = FALSE, bg = grey(0.6), fg = 0, add = TRUE, 
            xlim = c(m, M), ylim = c(0, 1))
    }
    
    stripchart(yF, pch = ">", method = "overplot", cex = 2, at = 0.55, axes = FALSE, xlim = c(m, M), ylim = c(0, 1), add = TRUE)
    mtext("IP", cex = 1.2, side = 2, las = 2, at = 0.5)
    stripchart(yR, pch = "<", method = "overplot", cex = 2, at = 0.45, col = 2, add = TRUE)
    
    abline(h = 0.35, lty = 3)
    if (addSe) {
        points(x@estimates$mu, rep(0.35, K), pch = "+", cex = 2)
        if (any(x@estimates$seMu != 0)) {
            points(x@estimates$mu - 2 * x@estimates$seMu, rep(0.35, K), pch = "[", cex = 1)
            points(x@estimates$mu + 2 * x@estimates$seMu, rep(0.35, K), pch = "]", cex = 1)
            segments(x@estimates$mu - 2 * x@estimates$seMu, rep(0.35, K), x@estimates$mu + 2 * x@estimates$seMu, rep(0.35, K), lwd = 1, lty = rep(1, K))
        }
    }
    
    if (length(cF) > 0) {
        stripchart(cF, pch = ">", method = "overplot", at = 0.25, cex = 2, add = TRUE, xlim = c(m, M), ylab = "Cont.", axes = FALSE)
    }
    if (length(cR) > 0) {
        stripchart(cR, pch = "<", method = "overplot", at = 0.15, cex = 2, col = 2, add = TRUE)
    }
    mtext("Cont.", cex = 1.2, side = 2, las = 2, at = 0.2)
    
    if (addNucleosome) {
        plot(c(m, M), c(0, 1), axes = FALSE, col = 0, ylim = c(0, 1), xlim = c(m, M), ylab = "")
        if (addSe) {
            symbols(x@estimates$mu, rep(0.5, K), rectangles = matrix(rep(c(147, 0.8), K), ncol = 2, byrow = TRUE), inches = FALSE, bg = grey(0.5 * pmin(se(x)/50, 
                1)), fg = 0, add = TRUE, xlim = c(m, M), ylim = c(0, 1))
        } else {
            symbols(x@estimates$mu, rep(0.5, K), rectangles = matrix(rep(c(147, 0.8), K), ncol = 2, byrow = TRUE), inches = FALSE, fg = 0, bg = 1, add = TRUE, xlim = c(m, 
                M), ylim = c(0, 1))
        }
        mtext("Nucl.", cex = 1.2, side = 2, las = 2, at = 0.5)
    }
})

#' @describeIn plot-pics Plot method for \code{picsError} and \code{segReads}
setMethod("plot", signature("picsError", "segReads"), function(x, y, addKernel = FALSE, main = NULL, ...) {
    par(oma = c(2.5, 5, 5, 5), mar = c(0, 5, 0, 0))
    layout(matrix(1:2, ncol = 1), heights = c(0.2, 0.1))
    
    yF <- y@yF
    yR <- y@yR
    cF <- y@cF
    cR <- y@cR
    map <- y@map
    m <- min(yF[1], yR[1]) - 100
    M <- max(tail(yF, 1), tail(yR, 1)) + 100
    
    stripchart(yF, pch = ">", method = "overplot", cex = 2, at = 0.5, add = FALSE, axes = FALSE, xlim = c(m, M), ylab = "Cont | Inp.", ylim = c(0, 1))
    stripchart(yR, pch = "<", method = "overplot", cex = 2, at = 0.5, col = 2, add = TRUE)
    
    abline(h = 0.35, lty = 3)
    
    # Add kernel density estimate
    if ((addKernel == TRUE) & (length(yF) > 1 & length(yR) > 1)) {
        dkF <- density(yF, bw = 75)
        dkR <- density(yR, bw = 75)
        plot(dkF, lty = 3)
        lines(dkR, col = 2, lty = 3)
    }
    
    if (length(cF) > 0) {
        stripchart(cF, pch = ">", method = "overplot", at = 0.2, cex = 2, add = TRUE, xlim = c(m, M), ylab = "Cont.", axes = FALSE)
    }
    if (length(cR) > 0) {
        stripchart(cR, pch = "<", method = "overplot", at = 0.2, cex = 2, col = 2, add = TRUE)
    }
})

#' @describeIn plot-pics Plot method for \code{picsList} and \code{segReadsList}
setMethod("plot", signature("picsList", "segReadsList"), function(x, y, regionIndex = NULL, addKernel = FALSE, addNucleosome = FALSE, addSe = TRUE, main = NULL, 
    ...) {
    if (is.null(main)) {
        setMain <- TRUE
    }
    if (is.null(regionIndex)) {
        regionIndex <- 1:length(x@List)
    }
    for (i in regionIndex) {
        if (setMain) {
            main <- paste(as.character(i), " (", y@List[[i]]@chr, ")", sep = "")
        }
        if (class(x@List[[i]]) != "picsError") {
            plot(x@List[[i]], y@List[[i]], addKernel = addKernel, addNucleosome = addNucleosome, addSe = addSe, main = main, ...)
        } else {
            plot(x@List[[i]], y@List[[i]], addKernel = addKernel, main = paste(as.character(i), " (", y@List[[i]]@chr, ")", sep = ""), ...)
            warning("Object of class picsError, no PICS density displayed")
        }
    }
})

#' Plot PICS FDR
#' 
#' This method plots a curve showing the FDR as a function of the PICS scores.
#' 
#' @param filter Filter to be passed to \code{picsFDR}.
#' @param h Argument passed to \code{abline}.
#' 
#' @param x A \code{picsList} object as returned by the function \code{PICS} run on the treatment data.
#' @param y A \code{picsList} object as returned by the function \code{PICS} run on the control data.
#' @param filter A list of ranges for filtering regions based on \code{PICS} parameters. By default filter is set to 'NULL' and all regions are used.
#' @param h A value between 0 and 1, representing the desired FDR. This simply draws a horizontal line at the given value.
#' @param ... Further graphical parameters passed to the generic function \code{plot}.
#' 
#' @seealso \code{PICS}
#' @name plot-FDR
#' @aliases plot,picsList,picsList-method
#' @export
setMethod("plot", signature("picsList", "picsList"), function(x, y, filter = NULL, h = 0.1, ...) {
    FDR <- picsFDR(x, y, filter = filter)
    arg <- list(...)
    par(mar = c(4, 4, 4.5, 4) + 0.1)
    if (length(arg$xlim) != 2) {
        xlim <- range(FDR[, 2])
    } else {
        xlim <- c(max(arg$xlim[1], min(FDR[, 2])), min(arg$xlim[2], max(FDR[, 2])))
    }
    plot(FDR[, 2], FDR[, 1], xlab = "score", ylab = "FDR", panel.first = grid(nx = 50), ...)
    xx <- FDR[FDR[, 2] > xlim[1] & FDR[, 2] < xlim[2], 2]
    yy <- FDR[FDR[, 2] > xlim[1] & FDR[, 2] < xlim[2], 3]
    xx <- xx[seq(1, length(xx), length.out = 10)]
    yy <- yy[seq(1, length(yy), length.out = 10)]
    axis(3, at = xx, labels = yy)
    mtext("# regions", side = 3, line = 3)
    FDRex <- FDR[FDR[, 1] > 0, ]
    notDup <- rev(!duplicated(rev(FDRex[, 1])))
    lines(FDRex[notDup, 2], FDRex[notDup, 1], col = 2, lty = 2, lwd = 1.5)
    abline(h = h, lw = 1.5, col = "grey")
})
