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
 \item{map}{\code{signature(x = ``pics'')}: subset PICS object.}
}
}


\author{
Charles Cheung, \email{cykc@interchange.ubc.ca} and
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