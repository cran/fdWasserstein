\name{Phoneme}
\alias{Phoneme}
\alias{logPeriodogram}
\docType{data}
\title{Phoneme data}
\description{
The dataset comprises 4509 log-periodograms computed from digitalized
speech frames. Each log-periodograms is of length 256, and is based on the
pronunciation of one of the following five phonemes:
"sh", "dcl", "iy", "aa" and "ao".   
}
\usage{
data(phoneme)
}
\format{
    \itemize{
        \item \code{logPeriodogram}: a 4509x256 matrix containing the log-periodograms.
        \item \code{Phoneme}: a vector of length 4509 containing the phonemes.
    }
}
\source{
The data set was downloaded from the "\emph{Elements of statistical
    learning}" website at \url{https://hastie.su.domains/ElemStatLearn/}
}
\references{
T. Hastie and R. Tibshirani and J. Friedman (2009) \emph{The elements of statistical learning: Data mining, inference and prediction}, 2nd edn, New York: Springer.
}
\examples{
data(phoneme)
old <- par(mfrow=c(3,2))
for (i in unique(Phoneme)) 
    matplot(t(logPeriodogram[Phoneme==i,]), type="l", 
            xlab="", ylab="", ylim=c(0,30), main=i)
par(old)
}
\keyword{datasets}

