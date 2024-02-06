## Partial spectral decomposition
## (only the eigenvalues, and corresponding eigenvectors,
##  greater than ..TOL.. are kept)
..TOL.. <- sqrt(.Machine$double.eps)

.eigen <- function(A) {
  e <- eigen(A, symmetric = TRUE)
  r <- which(e$values > ..TOL..)
  e$values <- e$values[r]
  e$vectors <- e$vectors[, r, drop = FALSE]
  e
}

## A^1/2
sqm <- function(A) {
  A <- .eigen(A)
  A$vectors %*% (sqrt(A$values) * t(A$vectors))
}

## A^-1/2
isqm <- function(A) {
  A <- .eigen(A)
  A$vectors %*% (1 / sqrt(A$values) * t(A$vectors))
}

## Wasserstein distance
dwstein <- function(A, B, B0.5) {
  if (missing(B0.5)) B0.5 <- sqm(B)
  max(
    0,
    sum(diag(A) + diag(B) - 2 *
          sqrt(pmax(
            0,
            eigen(B0.5 %*% A %*% B0.5, symmetric = TRUE, only.values = TRUE)$values
          )))
  )
}

## Wasserstein distance
dwasserstein <- function(A, B) dwstein(A, B, sqm(B))


## optimal map from A to B

## Version a) B^0.5 is available
.optGaussMap <- function(A, B0.5) B0.5 %*% isqm(B0.5 %*% A %*% B0.5) %*% B0.5

## Version b) Standard arguments
optGaussMap <- function(A, B) .optGaussMap(A, sqm(B))

## Frechet mean (Yoav algorithm)
## If max.iter=0 and gamma is missing,
## the Frechet mean with respect to the
## square root distance is returned
gaussBary <- function(sigma, w = rep(1, dim(sigma)[3]), gamma, sigma0.5,
                      max.iter = 30, eps = 1.0e-08, silent = max.iter == 0) {
  d <- dim(sigma)
  if (missing(sigma0.5)) sigma0.5 <- array(apply(sigma, 3, sqm), d)
  w <- rep(w / sum(w), rep(d[1] * d[2], d[3]))
  if (missing(gamma)) gamma <- crossprod(apply(w * sigma0.5, 1:2, sum))
  iter <- 0
  while (iter < max.iter) {
    T <- apply(sigma0.5, 3, function(S) .optGaussMap(gamma, S))
    Tm <- matrix(rowSums(w * T), d[1], d[1])
    g.old <- gamma
    gamma <- Tm %*% gamma %*% Tm
    ## cat(iter,sum(diag(Tm)),range(gamma-g.old),"\n")
    if (max(abs(gamma - g.old)) < eps) break
    iter <- iter + 1
  }
  if (!silent && (iter == max.iter)) warning("maximal number of iteration reached")
  list(gamma = gamma, iter = iter)
}
