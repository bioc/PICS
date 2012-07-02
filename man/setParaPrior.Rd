\name{setParaPrior}
\alias{setParaPrior}

\title{
	Function that returns a list of parameters that can be used as an argument of PICS and PING.
}

\description{
	This function takes from 0 to 6 parameters as argument, check if they are valid and returns a list to be used in a call to PICS or PING.
}

\usage{
	setParaPrior(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
}

\arguments{
	\item{xi}{ An \code{integer}, default=200. The average DNA fragment size.}
	\item{rho}{ An \code{integer}, default=1. A variance parameter for the average DNA fragment size distribution.}
	\item{alpha}{ An \code{integer}, default=20. First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model}
	\item{beta}{ An \code{integer}, default=40000. Second hyperparameter of the inverse Gamma distribution for sigma^2 in the PING model}
	\item{lambda}{ An \code{integer}, default=0. The precision of the prior for mu used for histone data.}
	\item{dMu}{ An \code{integer}, default=0. Our best guess for the distance between two neighboring nucleosomes.}
	\item{dataType}{ A \code{character}, The type of data you are  processing: default="TF" for transcription factor.}
}

\value{
	Returns a list of parameters to be used in PICS or PING.
}

\author{
Renan Sauteraud
}

\seealso{
	\code{PICS}
	\code{PING}
}

\keyword{functions}
