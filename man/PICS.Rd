\name{pics}
\alias{pics}
\alias{PICS}
\alias{pics-class}
\alias{as.list,AlignedRead-method}
\alias{code}
\alias{code,pics-method}
\alias{code,picsError-method}
\alias{code,picsList-method}
\alias{coerce}
\alias{coerce,picsList,RangedData-method}
\alias{coerce,picsList,data.frame-method}
\alias{coerce,AlignedRead,GenomeData-method}
\alias{coerce,RangedData,GenomeData-method}
\alias{coerce,data.frame,GenomeData-method}
\alias{coerce,data.frame,picsList-method}
\alias{coerce,data.frame,picsList-method}
\alias{density}
\alias{density,pics-method}
\alias{density,picsError-method}
\alias{score}
\alias{scoreForward}
\alias{scoreForward,pics-method}
\alias{scoreForward,picsList-method}
\alias{scoreForward,picsError-method}
\alias{scoreReverse}
\alias{scoreReverse,pics-method}
\alias{scoreReverse,picsList-method}
\alias{scoreReverse,picsError-method}
\alias{sigmaSqR}
\alias{sigmaSqR,pics-method}
\alias{sigmaSqR,picsList-method}
\alias{sigmaSqR,picsError-method}
\alias{sigmaSqF}
\alias{sigmaSqF,pics-method}
\alias{sigmaSqF,picsList-method}
\alias{sigmaSqF,picsError-method}
\alias{K}
\alias{K,pics-method}
\alias{K,picsList-method}
\alias{K,picsError-method}
\alias{plot,pics,segReads-method}
\alias{plot,picsError,segReads-method}
\alias{plot,picsList,segReadsList-method}
\alias{maxRange}
\alias{maxRange,pics-method}
\alias{maxRange,picsList-method}
\alias{maxRange,picsList-method}
\alias{maxRange,picsError-method}
\alias{minRange}
\alias{minRange,pics-method}
\alias{minRange,picsList-method}
\alias{minRange,picsList-method}
\alias{minRange,picsError-method}

% I think a few methods are missing, e.g. [[<- and [<-

\title{Estimation of binding site positions} 

\description{
This object contains Estimation of binding site positions and has the following  slots: 
segReadsList, dataType.
}

\usage{
PICS<-pics(segReadsList,dataType="TF")
}

\arguments{
\item{segReadsList}{This object contains segmentation of Genome}
\item{dataType}{ The type of data you are  processing: specified `TF' for transcription factor.}
}

\section{Methods}{
  \describe{
 \item{code}{\code{signature(x = ``pics'')}: return the error code for each list element (i.e. candidate region) of a PICS object. If the string is empty, there were no errors.}
 \item{plot}{\code{signature(x = ``pics'')}: Plot all regions in the PICS object. This might be long, and should only be used to plot a few regions, so subset the object before plotting.}
 \item{sigmaSqR}{\code{signature(x = ``pics'')}: return the variance parameter of the reverse (R) distribution for each binding event.}
 \item{sigmaSqF}{\code{signature(x = ``pics'')}: return the variance parameter of the forward (F) distribution for each binding event.}
 \item{score}{\code{signature(x = ``pics'')}: return the score for each binding event.} 
 \item{scoreF}{\code{signature(x = ``pics'')}: return the score  of the forward (F)  for each binding event.}
 \item{scoreR}{\code{signature(x = ``pics'')}: return the score  of the forward (R)  for each binding event.}
\item{maxRange}{\code{signature(x = ``pics'')}: return the range maximum.}
\item{minRange}{\code{signature(x = ``pics'')}: return the range minimal.}
 \item{K}{\code{signature(x = ``pics'')}: subset PICS object.}
 \item{density}{\code{signature(x = ``pics'')}: return the density for each binding event.}
}
}

\author{
  Xuekui Zhang, \email{ubcxzhang@gmail.com} and
  Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
  Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\seealso{
\code{\link{pics}}
}
\keyword{models}
