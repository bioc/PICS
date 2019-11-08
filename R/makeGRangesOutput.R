#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges strand seqnames
#' @importFrom IRanges IRanges coverage
#' 
#' @export
makeGRangesOutput <- function(obj, type = "fixed", length = 100,
                              filter = list(delta = c(0, Inf), se = c(0, Inf), 
                                            sigmaSqF = c(0, Inf), sigmaSqR = c(0, Inf), score = c(0, Inf))){
    nSe <- 3
    if (type != "wig") {
        mu <- mu(obj)
        delta <- delta(obj)
        se <- se(obj)
        seF <- seF(obj)
        seR <- seR(obj)
        sf <- sigmaSqF(obj)
        sr <- sigmaSqR(obj)
        chromosome <- chromosome(obj)
        score <- score(obj)
        
        ## Filter regions with small deltas over a region
        if (!is.null(filter)) {
            ### Filter based on delta
            ind1 <- delta > filter$delta[1] & delta < filter$delta[2]
            ind2 <- sf > filter$sigmaSqF[1] & sf < filter$sigmaSqF[2]
            ind3 <- sr > filter$sigmaSqR[1] & sr < filter$sigmaSqR[2]
            ind5 <- is.finite(score) & score > filter$score[1] & score < filter$score[2]
            ind <- ind1 & ind2 & ind3 & ind5
            
            if (!is.null(filter$se)) {
                ind4 <- se > filter$se[1] & se < filter$se[2] & seF > filter$se[1] & seF < filter$se[2] & seR > filter$se[1] & seR < filter$se[2]
                ind <- ind & ind4
            }
        } else {
            ind <- is.finite(score)
        }
    }
    if (type == "bed") {
        score <- (score(obj))[ind]
        ord <- order(-score)
        # Order the score
        score <- score[ord]
        start <- (mu - delta/2 - nSe * seF)[ind]
        end <- (mu + delta/2 + nSe * seR)[ind]
        start <- start[ord]
        end <- end[ord]
        # chrom<-(paste('chr', chromosome, sep=''))[ord]
        chrom <- (paste("", chromosome, sep = ""))[ind]
        chrom <- chrom[ord]
    } else if (type == "ci") {
        score <- (score(obj))[ind]
        ord <- order(-score)
        score <- score[ord]
        start <- (mu - nSe * se)[ind]
        end <- (mu + nSe * se)[ind]
        start <- start[ord]
        end <- end[ord]
        chrom <- (paste("", chromosome, sep = ""))[ind]
        chrom <- chrom[ord]
    } else if (type == "fixed") {
        score <- (score(obj))[ind]
        ord <- order(-score)
        score <- score[ord]
        start <- (mu - length)[ind]
        end <- (mu + length)[ind]
        start <- start[ord]
        end <- end[ord]
        chrom <- (paste("", chromosome, sep = ""))[ind]
        chrom <- chrom[ord]
        strand <- NULL
    } else if (type == "wig") {
        temp <- wigDensity(obj, strand = "*", step = 10, sum = TRUE, filter = filter, scale = TRUE)
        chrom <- temp$chr
        start <- temp$x
        end <- temp$x + 9
        score <- temp$density
        strand = "*"
    }
    
    if (length(score) == 0) {
        gr = NULL
    } else {
        ranges <- IRanges(as.integer(start), as.integer(end))
        if (type == "bed" | type == "ci" | type == "fixed") {
            names(ranges) <- paste("pics", 1:(length(score)), sep = "")
            gr <- GRanges(seqnames = chrom, ranges = ranges, strand = '*', score)
        } else {
            gr <- GRanges(seqnames = chrom, ranges = ranges, strand = strand, score)
            print("Removing overlapping binding events")
            overlap <- as.matrix(findOverlaps(ranges(gr)))
            overlap <- split(overlap[, 1], overlap[, 2])
            scL <- score(gr)
            maxScores <- lapply(overlap, function(x) {
                max(scL[x])
            })
            maxScscores <- as.numeric(unlist(maxScores))
            idx <- which(maxScores == scL)
            gr <- gr[idx, ]
        }
    }
    return(gr)
}
