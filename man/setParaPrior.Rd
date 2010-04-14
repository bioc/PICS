\name{setParaPrior} 
\alias{setParaPrior} 
\alias{paraPriorTF}
\alias{paraPriorH}

\title{Set convergence parameters of the EM algorithm} 

\description{ 
This function can be used to change the internal PICS parameters for the prior distribution. This function should only be called if you really now what you are doing! In particular, you may want to specify the average DNA fragment size for your sample by changing the `xi' parameter.
}
\usage{
setParaPrior(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")

}
\arguments{
  \item{xi}{Our best guest for the average DNA fragment size.}
  \item{rho}{A variance parameter for the average DNA fragment size distribution.}
  \item{alpha}{First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model.}
  \item{beta}{First hyperparameter of the inverse Gamma distribution for sigma^2 in the PICS model.}
  \item{lambda}{The precision of the prior for mu used for histone data.}
  \item{dMu}{Our best guess for the distance between two neighboring nucleosomes.}
  \item{dataType}{A character string equal to either `H' or `TF', `H' for histone and `TF' for transcription factors. `H' is not yet supported}
}
%\details{
%}
\value{
No value returned. The function simply modifies the internal variables `paraPriorTF' if dataType=`TF' and `paraPriorH' if dataType='H'.
}
\references{
X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S. Jones, and R. Gottardo, ``PICS: Probabilistic Inference for ChIP-seq'' arXiv, 0903.3206, 2009. To appear in Biometrics.
}
\author{Xuekui Zhang and Raphael Gottardo <\email{raphael.gottardo@ircm.qc.ca}>}

\seealso{setParaEM}
\examples{
# Using xi=200 for the average DNA fragment size
setParaPrior(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
# Using xi=150 for the average DNA fragment size
setParaPrior(xi=150,rho=1,alpha=20,beta=40000,lambda=0,dMu=200,dataType="TF")
}

\keyword{data}