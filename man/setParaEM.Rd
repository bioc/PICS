\name{setParaEM} 
\alias{setParaEM} 
\alias{paraEMH} 
\alias{paraEMTF} 
\title{Set convergence parameters of the EM algorithm} 

\description{ 
This function can be used to change the internal PICS parameters for the EM algorithm. This function should only be called if you really now what you are doing!.}

\usage{
setParaEM(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType="TF")
}

\arguments{
  \item{minK}{An integer value. The minimum number of binding events per region.}
  \item{maxK}{An integer value. The maximum number of binding events per region.}
  \item{tol}{The tolerance for the EM algorithm}
  \item{B}{An integer value. The maximum number of iterations to be used.}
  \item{mSelect}{A character string specifying the information criteria to be used when selecting the number of binding events.}
  \item{mergePeaks}{A logical value stating whether overlapping binding events should be picked.}
  \item{mapCorrect}{Should mappability profiles be incorporated in the estimation, that is missing reads estimated.}
  \item{dataType}{A character string equal to either 'H' or 'TF', 'H' for histone and 'TF' for transcription factors.}
}


\value{
No value returned. The function simply modifies the internal variables 'paraEMTF' if dataType='TF' and 'paraEMH' if dataType='H'.
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, "PICS: Probabilistic Inference for ChIP-seq," arXiv, 0903.3206, 2009.
}
\author{Xuekui Zhang and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>}

\seealso{setParaPrior}
\examples{
# Using mSelect="BIC"
setParaEM(minK=1,maxK=15,tol=1e-4,B=100,mSelect="BIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType="TF")
# Using mSelect="AIC"
setParaEM(minK=1,maxK=15,tol=1e-4,B=100,mSelect="AIC",mergePeaks=TRUE,mapCorrect=TRUE,dataType="TF")
}
\keyword{data}