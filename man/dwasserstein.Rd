\name{dwasserstein}
\alias{dwasserstein}

\title{
2-Wasserstein distance
}
\description{
  Computes the 2-Wasserstein distance between the (covariance) matrices A and B.
}
\usage{
  dwasserstein(A, B)
}

\arguments{
  \item{A,B}{
  Two symmetric positive semi-definite matrices.
}
}

\value{
A numeric object with the 2-Wasserstein distance of A and B. 
}
\references{
Masarotto, V., Panaretos, V.M. & Zemel, Y. (2019) "Procrustes Metrics on
Covariance Operators and Optimal Transportation of Gaussian Processes",
\emph{Sankhya A} \strong{81}, 172-213 \doi{10.1007/s13171-018-0130-1}
}
\author{
Valentina Masarotto, Guido Masarotto} 

\seealso{
\code{\link{gaussBary}}
}
\examples{
n <- 10
matrices <- rWishart(2,n,diag(n))
A <- matrices[,,2]
B <- matrices[,,1]
dwasserstein(A,B)
dwasserstein(A, 10*crossprod(B))
}
