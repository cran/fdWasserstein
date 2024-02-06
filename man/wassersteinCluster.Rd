\name{wassersteinCluster}
\alias{wassersteinCluster}
\alias{trimmedAverageSilhouette}

\title{
Soft clustering of covariance operators.
}
\description{
%wassersteinCluster` is the main function. 
Computes the soft cluster solutions for different values of the number of clusters K. 
}
\usage{
wassersteinCluster(data, grp, 
                   kmin = 2, kmax = 10, 
                   E = -0.75 * (0.95 * log(0.95) + 
                        0.05 * log(0.05)) + 0.25 * log(2), 
                   nstart = 5, nrefine = 5, ntry = 0, 
                   max.iter = 20, tol = 0.001, 
                   nreduced = length(unique(grp)), 
                   nperm = 0, 
                   add.sigma = FALSE, 
                   use.future = FALSE, verbose = TRUE)

trimmedAverageSilhouette(a, plot = TRUE)
}

\arguments{
  \item{data}{
A N times M matrix containing the N sample curves; M denotes the number of points of the grid on which the curves are available.
}
  \item{grp}{
A vector or factor of length N; a covariance operator is estimated for
each level of grp.
}
  \item{kmin,kmax}{
A pair of integer defining the desired number of clusters. A solution is computed for K=kmin,...,kmax.
}
  \item{E}{
The desired average entropy.
}
  \item{nstart, nrefine, ntry}{
The integers used during the initialization search. If ntry=0, then `ntry` is set to `round(1+N/K)`.
}
  \item{max.iter}{
Maximum number of block descend iterations.
}
  \item{tol}{
Iterations stop when the relative decrease of the objective function in two consecutive iterations is less than `tol`.
}
  \item{nreduced}{
The number of covariances used to estimate the cluster barycenters.
}
  \item{nperm}{
The number of permutation used to approximate the reference distribution of max TASW.
}
  \item{add.sigma}{
Should the sample covariances be returned?
}
  \item{use.future}{
Use or not use package 'future' to parallelize the computation? See note.
}
  \item{verbose}{
If `verbose==TRUE`, information on the progress of the optimization are shown.
}
  \item{a}{
A list returned by `wassersteinCluster`.
}
  \item{plot}{
If `plot==TRUE`, the TASW profile is plotted.
}
}


\value{
`wassersteinCluster` returns a list of length kmax-kmin+1. 
The ith element is a list describing the cluster solution obtained for k=kmin+i-1, and containing: 
\item{K,E,eta}{the number of cluster, the average entropy and the corresponding value of `eta`;}
\item{w}{the N times K soft partition matrix;}
\item{g}{a M times M times K array with the cluster barycenters;}
\item{d}{a N times K matrices containing the distances between the N sample covariances and the K cluster barycenters;}
\item{obj}{`obj`: the minimum value of the objective function.}
The list may have the following attributes:
\item{df}{the degree of freedom of the sample operators (a vector). Always present.}
\item{sample.covariances}{a list contaning the sample operators (as a
     3-dimensional array); only present if add.sigma=TRUE;}
\item{tasw.test}{a list containing the value of maxTASW computed from the data (a scalar), the nperm values of of maxTASW obtained by permutation (a vector), and the corresponding p-value (a scalar); only present if nperm>0.}

`trimmedAverageSilhouette` returns a numeric vector with the TASW values.

}
\details{
See Masarotto & Masarotto (2023) for the algorithm details.
}
\references{
Masarotto, V. & Masarotto, G. (2023) "Covariance-based soft clustering
of functional data based on the Wasserstein-Procrustes metric", 
\emph{Scandinavian Journal of Statistics}, \doi{10.1111/sjos.12692}.
}
\author{
Valentina Masarotto, Guido Masarotto
}
\note{
To distribute the computation on more than a cpu 
\enumerate{
\item install the package 'future'
\item  execute in the R session
\itemize{
\item library(future)
\item plan(multissession)
}
}
For more options, see the future's documentation
}

\examples{
\donttest{
# Example phoneme.R (simplified) from https://doi.org/10.1111/sjos.12692. 
data(phoneme)
# resampling the log-periodograms
# 15 sample covariances for each phoneme
set.seed(12345)
nsubsamples <- 15
n <- 40
gg <- unique(Phoneme)
nphonemes <- length(gg)
N <- n*nsubsamples*nphonemes
M <- NCOL(logPeriodogram)
X <- matrix(NA, N, M)
gr <- integer(N)
r <- 1
first <- 1
last <- n
for (l in gg) {
  for (i in 1:nsubsamples) {
    X[first:last, ] <- logPeriodogram[sample(which(Phoneme==l),n), ]
    gr[first:last] <- r
    r <- r+1
    first <- first+n
    last <- last+n
  }
}
# soft clustering
a <- wassersteinCluster(X, gr)
# how many cluster?
trimmedAverageSilhouette(a)
# the membership weigths show that the
# algorithm reconstructed the five phoneme
w <- ts(a[[4]]$w)
colnames(w) <- paste("Cluster", 1:5)
plot(w, xlab="Sample covariances", main="")
}
}

\keyword{multivariate}
\keyword{cluster}
\concept{functional analysis}






