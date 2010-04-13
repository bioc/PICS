\name{segReadsList}
\docType{class}
\alias{segReadsList}
\alias{segReadsList-class}
\alias{[,segReadsList,ANY,ANY-method}
\alias{[[,segReadsList,ANY,ANY-method}
\alias{length,segReadsList-method}

\title{Segment the genome into candidate regions} 

\description{
Pre-process bidirectional aligned reads data from a single ChIP-Seq experiment to detect candidate regions with a minimum number of forward and reverse reads. These candidate regions will then be processed by PICS.}

\usage{
segReadsList<-new(List,paraSW,N,Nc)
}

\arguments{
\item{List}{This object contains an ExpressionSet}
\item{paraSW}{ String containing the genome name used (vector).}
\item{N}{ String containing the name of chromosome used (vector).}
\item{Nc}{ String containing the Position of the sequences (vector).}
}


\section{Methods}{
  \describe{
 \item{[}{\code{signature(x = "pics")}: subset gadem object.}
\item{[}{\code{signature(x = "pics")}: subset gadem object.}
	}
}


\section{Methods}{
  \describe{
 \item{length}{\code{signature(x = "pics")}: subset PICS object.}
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