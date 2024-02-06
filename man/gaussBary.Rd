\name{gaussBary}
\alias{gaussBary}

\title{
Wasserstein barycenter between Gaussian Processes
}
\description{
Computes the Frechet mean between covariance operators with respect to the Procrustes metrics (equivalently, a Wasserstein barycenter of centered Gaussian processes with corresponding covariances) via steepest gradient descent.
}
\usage{
gaussBary(sigma, w = rep(1, dim(sigma)[3]), gamma, sigma0.5, 
  max.iter = 30, eps = 1e-08, silent = max.iter == 0)
}

\arguments{
  \item{sigma}{
An MxMxK array containing the K covariances. 
}
  \item{w}{
Optional. A vector of weights of length K. If missing, each matrix is given equal weight 1. 
}
  \item{gamma}{
Optional. Initialisation point for the gradient descent algorithm. 
}
  \item{sigma0.5}{
    Optional. An array containing the square roots of the matrices in
    sigma if available. The square roots are computed by
    \code{gaussBary} if \code{sigma0.5} is missing.
}
  \item{max.iter}{
Maximum number of gradient descent iterations.
}
  \item{eps}{
Iterations stop when the relative decrease of the objective function in two consecutive iterations is less than `eps`. 
}
  \item{silent}{
If \code{FALSE} returns a warning if maximal number of iteration is reached. 
}
}


\value{
A list of 2 containing:
\item{gamma}{The MxM Frechet mean.}
\item{iter}{Number of iterations needed to reach convergence, numeric.}
}
\references{
Masarotto, V., Panaretos, V.M. & Zemel, Y. (2019) "Procrustes Metrics on
Covariance Operators and Optimal Transportation of Gaussian Processes",
\emph{Sankhya A} \strong{81}, 172-213 \doi{10.1007/s13171-018-0130-1}
}
\author{
Valentina Masarotto, Guido Masarotto
}
\note{
We thank Yoav Zemel for the first version of the code. 
}

\examples{
M <- 5
K <- 4

sigma <- rWishart(M, df = K, Sigma = diag(K))

gaussBary(sigma) 
}

\keyword{multivariate}
\concept{functional analysis}
