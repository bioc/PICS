\name{segReads}
\alias{segReads}
\alias{segReads-class}
\alias{map}
\alias{map,segReads-method}
\alias{map,segReadsList-method}
\alias{\%in\%,ANY,segReads-method}


\title{Segment the genome into candidate regions} 

\description{
Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PICS.
}

\usage{
segReads<-new(yF, yR, cF, cR, map, chr)
}

\arguments{
\item{yF}{This object contains an ExpressionSet}
\item{yR}{ String containing the genome name used (vector).}
\item{cF}{ String containing the name of chromosome used (vector).}
\item{cR}{ String containing the Position of the sequences (vector).}
\item{map}{ String containing the copy number of sequence (vector).}
\item{chr}{ String containing the expresion data of enriched region (matrix with n column.}

}

\section{Methods}{
  \describe{
 \item{map}{\code{signature(x = "pics")}: subset PICS object.}
}
}


\author{
Charles Cheung, \email{cykc@interchange.ubc.ca} and
Raphael Gottardo, \email{raphael.gottardo@ircm.qc.ca}
Arnaud Droit, \email{arnaud.droit@ircm.qc.ca}
}

\references{
W. E. Johnson, Li, W., Meyer, C. A., Gottardo, R., Carroll, J. S., Brown, M., and Liu, X. S. (2006). Model-based analysis of tiling-arrays for ChIP-chip. PNAS 103:12457-12462.
}

\seealso{
\code{\link{pics}}
}

\keyword{models}