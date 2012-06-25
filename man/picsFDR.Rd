\name{picsFDR}
\alias{picsFDR}


\title{Compute the global FDR}
\description{
  Calculate an estimate of the FDR for PICS. This calculation requires control data (e.g. from an input DNA sample). 
}
\usage{
picsFDR(picsIP,picsCont,filter=list(delta=c(0,Inf),se=c(0,Inf),sigmaSqF=c(0,Inf),sigmaSqR=c(0,Inf)))
}

\arguments{
  \item{picsIP}{An object of class 'picsList' as returned by 'PICS' when run on IP/Control data.}
  \item{picsCont}{An object of class 'picsList' as returned by 'PICS' when run on Control/IP data.}
  \item{filter}{A list of filters to be used before computing the FDR. By default all regions are included. See details for more info on how to specify the filters.}
  
}

\value{
A dataframe with three columns corresponding to the estimated FDR, the score, and the number of regions. 
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\author{ Xuekui Zhang, Arnaud Droit <\email{arnaud.droit@crchuq.ualaval.ca}> and Raphael Gottardo <\email{rgottard@fhcrc.org}>}
\seealso{
\code{\link{pics}}
}
\examples{
\dontrun{
  # Segment the reads in order to identify candidate regions
   segC<-segmentReads(RD, minReads=3,invert=TRUE)
  # Use the serial version of PICS
  picsC<-PICS(segC,dataType='TF')
	plot(pics,picsC,xlim=c(0,50),ylim=c(0,.2),filter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500)),type="l")  
}
}

\keyword{models}

