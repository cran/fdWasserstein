\name{tangentPCA}
\alias{tangentPCA}
\title{
Tangent space principal component analysis 
}
\description{
The function performs a standard PCA of K covariances after projecting them
on the tangent space at their Wasserstein barycenter. Rationale and
details are given in Masarotto, Panaretos & Zemel (2019, 2022)
}
\usage{
tangentPCA(sigma, max.iter=30)
}
\arguments{
  \item{sigma}{
An MxMxK array containing the K covariances. 
}
\item{max.iter}{
Maximum number of gradient descent iterations used to compute the
Wasserstein barycenter of the covariances in sigma.
}
}
\value{
A standard \code{\link{prcomp}} object with added a MxMxK array
containing the eigenvectors projected back to the covariances space.
}
\references{
Masarotto, V., Panaretos, V.M. & Zemel, Y. (2022) "Transportation-Based Functional ANOVA and PCA for Covariance Operators", \emph{arXiv},  \url{https://arxiv.org/abs/2212.04797}
}
\author{
Valentina Masarotto
}

\seealso{
\code{\link{gaussBary}}, \code{\link{prcomp}}
}
\examples{
\donttest{
# Example taken from https://arxiv.org/abs/2212.04797  . 
data(phoneme)
# resampling the log-periodograms
# 12 sample covariances for each phoneme
# each estimated on 50 curves 
set.seed(12345)
nsubsamples <- 12
n <- 50
gg <- unique(Phoneme)
nphonemes <- length(gg)
K <- n*nsubsamples*nphonemes
M <- NCOL(logPeriodogram)
Sigma <- array(dim=c(M, M, nphonemes*nsubsamples))
r <- 0
for (l in gg) {
  for (i in 1:nsubsamples) {
      r <- r+1
      Sigma[,,r] <- cov(logPeriodogram[sample(which(Phoneme==l),n), ])
  }
}
pca <- tangentPCA(Sigma, max.iter=3)
summary(pca)
plot(pca)
# See https://arxiv.org/abs/2212.04797 for the interpretation
# of the figure
pairs(pca$x[,1:5], col=rep(1:nphonemes, rep(nsubsamples, nphonemes)))
}
}

\keyword{multivariate}
\concept{functional analysis}
