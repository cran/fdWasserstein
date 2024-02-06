\name{wassersteinTest}
\alias{wassersteinTest}

\title{
A permutation or bootstrap test based on optimal transport maps. 
}
\description{
The main function performs a k-sample permutation- or bootstrap-based test to check the equality of covariance operators. More specifically, given a sample of N functional curves belonging to K different populations, each characterized by its own covariance operators, the test aims to check the null hypothesis \eqn{\Sigma_1 = \dots = \Sigma_K}{Sigma1 = .... = SigmaK} versus the alternative that at least one operator is different. The test leverages on the equivalence between covariance operators and centered Gaussian processes. In the default version, in order to test the null, the test builds optimal transport maps from the sample to the Wasserstein barycenter of the processes. Successively, it contrasts these maps to the identity operator, as explained in Masarotto, Panaretos & Zemel (2022). However, argument "statistics" allows to base the test directly on the Wasserstein distance between covariance operators, rather than on optimal maps.
}
\usage{
wassersteinTest(data, grp, B = 1000, 
                statistic = c("transport", "distance"), 
                type = c("permutation", "bootstrap"), 
                r = c("HS", "trace", "operator"),
                align = TRUE, 
                use.future = FALSE,
                iter.bary = 10)
}

\arguments{
  \item{data}{
A N times M matrix containing the N sample curves; M denotes the number of points of the grid on which the curves are available.
}
\item{grp}{
Labels that identify which population each curve belongs to.
}
  \item{B}{
Number of permutations or bootstrap replications. If missing, B=1000. 
}
  \item{statistic}{
Whether the test is based on the transport maps or directly on the Wasserstein distance. Default is transport.  
}
  \item{type}{
Whether the test is permutation or bootstrap based. 
}
  \item{r}{
  Which norm is used to contrast the test statistics to 0 (used only if statistics="transport"). If r="HS" the Hilbert-Schmidt norm is used, if r="trace" the trace (nuclear) norm is used, if r="operator", the operator norm is used. Default is r="HS".
  }
  \item{align}{
If 'align=TRUE', the curves are centered around their mean. Default is TRUE.
}
  \item{use.future}{
Use or not use package 'future' to parallelize the computation? See note.
}
  \item{iter.bary}{
After how many iterations the gradient descent algorithm to compute the barycenter stops. 
}
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
For more options, see the future's documentation.
}

\value{
A list of three returning:
\item{stat}{Observed value of the test statistics}
\item{p.value}{The p-value indicating the significance level of the test}
\item{trep}{Value of the test statistics for each of the B permutation}
}
\references{
Masarotto, V., Panaretos, V.M. & Zemel, Y. (2022) "Transportation-Based Functional ANOVA and PCA for Covariance Operators", \emph{arXiv},  \url{https://arxiv.org/abs/2212.04797}
}
\author{
Valentina Masarotto, Guido Masarotto
}

\examples{
n = 20 
size <- 10
covariances <- rWishart(2,size,diag(size))
A <- covariances[,,1]
B <- covariances[,,2]

# Two groups, each with one covariance. Creates n Gaussian data for each covariance.
# more generally, we could have two groups each with "g_i" covariances in them
g1 <- g2 <- 1
grp <- rep(1:(g1+g2),rep(n,g1+g2))

data <- rbind(matrix(rnorm(n*NCOL(A)),n*g1)\%*\%A,
             matrix(rnorm(n*NCOL(B)),n*g2)\%*\%B)
wassersteinTest(data,grp, B=100,r="HS")$p.value
\donttest{
data(phoneme)
wassersteinTest(logPeriodogram, Phoneme, B=100,r="HS")$p.value
}
}

\keyword{multivariate}
\keyword{htest}
\concept{functional analysis}

