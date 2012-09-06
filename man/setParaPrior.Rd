\name{setParaPrior}
\alias{setParaPrior}

\title{
	Function that returns a list of parameters that can be used as an argument of PICS.
}

\description{
	This function takes from 0 to 6 parameters as argument, check if they are valid and returns a list to be used in a call to PICS.
}

\usage{
	setParaPrior(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200, dataType=NULL)
}

\arguments{
	\item{xi}{ An \code{integer}, default=200. The average DNA fragment size.}
	\item{rho}{ An \code{integer}, default=1. A variance parameter for the average DNA fragment size distribution.}
	\item{alpha}{ An \code{integer}, default=20. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
	\item{beta}{ An \code{integer}, default=40000. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
	\item{lambda}{ An \code{integer}, default=0. The precision of the prior for mu used for histone data.}
	\item{dMu}{ An \code{integer}, default=0. Our best guess for the distance between two neighboring nucleosomes.}
	\item{dataType}{ A character string. If a valid dataType is specified, use our suggested parameters.}
}

\value{
	Returns a list of 6 parameters to be used in PICS.
}

\author{
Renan Sauteraud
}

\seealso{
	\code{PICS}
}

\examples{
	# set prior for PICS data
	paraPrior<-setParaPrior(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200)
	# set prior for chip-seq data using our selected default parameters
	paraPrior<-setParaPrior(dataType="chip-seq")
	
}

\keyword{functions}
