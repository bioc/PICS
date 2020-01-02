#' @name segmentPICS
#' @title Segment the genome into candidate regions
#' 
#' @description 
#' Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment
#'  to detect candidate regions with a minimum number of forward and reverse 
#'  reads. These candidate regions will then be processed by PICS.
#'
#' @param data A \code{GRanges} object containing the IP reads. See details for 
#' more information on how to set up the data.
#' @param dataC A \code{GRanges} object containing the control reads. Set to 
#' NULL by default, i.e. no control.
#' @param map A \code{GRanges} object containing the mappability profiles. Set 
#' to NULL by default, i.e. no profiles.
#' @param minReads A \code{numeric}. The minimum number of F/R reads to be 
#' present in the sliding window.
#' @param minReadsInRegion A \code{numeric}. The minimum number of F/R reads to
#'  be present in the region.
#' @param jitter	A \code{logical} value stating whether some noise should be 
#' added to the read locations. This is recommended if the read positions have lots of duplicates.
#' @param dataType A \code{character}. Type of experiment. "TF" or "H".
#' @param maxLregion A \code{numeric}. The maximum length.
#' @param minLregion A \code{numeric}. The minimum length.
#'
#' @return An object of class \code{segReadsList} containing the results for all
#' pre-processed regions.
#' 
#' @examples 
#' # Read data
#' path<-system.file("extdata",package="PICS")
#' ## Note that the col name for the chromosome needs to be space and not chr
#' dataIP <- read.table(file.path(path, "Treatment_tags_chr21_sort.bed"), header=TRUE,
#'                      colClasses = c("factor","integer","integer","factor"))
#' dataIP <- as(dataIP, "GRanges")
#' 
#' dataCont <- read.table(file.path(path, "Input_tags_chr21_sort.bed"), header=TRUE,
#'                        colClasses = c("factor","integer","integer","factor"))
#' dataCont <- as(dataCont, "GRanges")
#' 
#' map <- read.table(file.path(path, "mapProfileShort"), header=TRUE,
#'                   colClasses = c("factor","integer","integer","NULL"))
#' map <- as(map, "GRanges")
#' seg <- segmentPICS(dataIP, dataC = dataCont, map = map, minReads = 1)
#'
#' @references 
#' X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo,
#' “PICS: Probabilistic Inference for ChIP-seq” arXiv, 0903.3206, 2009.
#'
#' @seealso segReadsList
#'
#' @importFrom GenomicAlignments readGAlignments
#' @export
segmentPICS<-function(data, dataC=NULL, map=NULL, minReads=2, minReadsInRegion=3,
                      jitter=FALSE, dataType="TF",maxLregion=0,minLregion=100){
  #Paras depends on the datatype
  step=20
  if(dataType=="TF") width=250
  if(dataType=="H")  width=150
  
  newSet<-segReadsGeneric(data, dataC=dataC, map=map, minReads=minReads, minReadsInRegion=minReadsInRegion, jitter=jitter, maxLregion=maxLregion,minLregion=minLregion, step=step, width=width, package="PICS")
  return(newSet)
}

#' Pre-process bam files
#' 
#' @description 
#' Reads a bam file using \code{Rsamtools} and extract the reads for each chromosome.
#' 
#' @param bamFile A \code{character}. The name of the .bam file to read.
#' @param chr A \code{character}. An optional character string. If specified, 
#'  only the selected chromosome will be returned. Speed up the computation.
#' @param PE A \code{logical}. Set to \code{TRUE} for paired-end sequencing
#'  data.
#' @param verbose A \code{logical}. Print additional information about the data.
#' 
#' @return A \code{GRanges} of all the reads for each chromosome.
#' 
#' @note
#' The user might encounter a memory allocation error when using bam 
#'  files of bigger sizes. Splitting the file by chromosome before calling \code{bam2gr}
#'  will solve this issue.
#'  
#' For Paired-End data, non matched reads are discarded.
#' 
#' @seealso segmentPICS
#' 
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicRanges elementMetadata
#' @export
bam2gr<-function(bamFile, chr=NULL, PE=FALSE, verbose=FALSE){
  paras <- ScanBamParam(what=c("qname", "rname", "strand","mapq"), flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE))
  bga<-readGAlignments(bamFile, use.names=FALSE, param=paras)
  if(verbose){
    cat(length(bga)," Reads in '",bamFile,"'","\n", sep="")
  }
  hiQScoreIdx<-which(elementMetadata(bga)$mapq>10)
  if(verbose){
    cat(length(bga)-length(hiQScoreIdx)," Reads with low quality scores filtered out","\n")
  }
  bga<-bga[hiQScoreIdx]#filter out elt with low quality scores

  if(isTRUE(PE)){
    qname<-elementMetadata(bga)$qname
    qname<-substring(qname,15)
    qname<-gsub(pattern="/3", replacement="", qname)
    elementMetadata(bga)$qname<-qname
    #merge pairs
    asdf<-as(bga, "data.frame")
    if(verbose){
      df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide")
    } else{
      suppressWarnings(df<-reshape(asdf, timevar="strand", idvar="qname", direction="wide"))
    }   
    df2<-df[,c("start.+", "end.-","rname.+")]
    colnames(df2)<-c("start", "end", "chr")
    rownames(df2)<-df[,"qname"]
    badReads<-which(df2$start>df2$end)
    if(length(badReads)>0){
      df2<-df2[-badReads,]
    }   
    #Split PE and SE
    idx <- is.na(df2[,c("start","end")])
    reads <- df2[!(idx[,1]|idx[,2]),]
    yFm <- df2[idx[,2],] #Forward reads, missing matching reverse
    yRm <- df2[idx[,1],] #Reverse reads
    gr<-GRanges(ranges=IRanges(start=reads$start, end=reads$end), strand="*", seqnames=reads$chr)
  } else{
    gr<-GRanges(ranges=IRanges(start=start(bga), end=end(bga)), strand=strand(bga), seqnames=seqnames(bga))
  }
  return(gr)
}

