\name{fdWasserstein-package}
\alias{fdWasserstein-package}
\alias{fdWasserstein}
\docType{package}
\title{
\packageTitle{fdWasserstein}
}
\description{
A package containing functions developed to support statistical analysis on functional covariance operators. 
In particular, 
\itemize{
    \item Function \code{\link{dwasserstein}} computes the
    Wasserstein-Procrustes distance between two covariances.

    \item Function \code{\link{gaussBary}} computes the Frechet mean of
    K covariances with respect to the Procrustes metrics
    (equivalently, the Wasserstein barycenter of centered Gaussian
    processes with corresponding covariances) via steepest gradient
    descent. See Masarotto, Panaretos & Zemel (2019).

    \item Function \code{\link{tangentPCA}} performs the tangent space
    principal component analysis considered in Masarotto, Panaretos &
    Zemel (2022). 

    \item Function \code{\link{wassersteinTest}} lets to test the null
    hypothesis that K covariances are equal using the methodology suggested by
    Masarotto, Panaretos & Zemel (2022).

   \item Function \code{\link{wassersteinCluster}} implements the soft
   partion procedure proposed by Masarotto & Masarotto (2023).
}
}
\author{
\packageAuthor{fdWasserstein}

Maintainer: \packageMaintainer{fdWasserstein}
}
\references{
Masarotto, V., Panaretos, V.M. & Zemel, Y. (2019) "Procrustes Metrics on
Covariance Operators and Optimal Transportation of Gaussian Processes",
\emph{Sankhya A} \strong{81}, 172-213 \doi{10.1007/s13171-018-0130-1}

Masarotto, V., Panaretos, V.M. & Zemel, Y. (2022) "Transportation-Based Functional ANOVA and PCA for Covariance Operators", \emph{arXiv},  \url{https://arxiv.org/abs/2212.04797}

Masarotto, V. & Masarotto, G. (2023) "Covariance-based soft clustering
of functional data based on the Wasserstein-Procrustes metric", 
\emph{Scandinavian Journal of Statistics}, \doi{10.1111/sjos.12692}.
}

