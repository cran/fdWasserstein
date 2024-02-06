tangentPCA <- function(sigma, max.iter = 30) {
    size <- dim(sigma)
    bary <- gaussBary(sigma, max.iter = max.iter, silent = TRUE)$gamma
    mydat.tg <- array(dim = size)

    # lifting to tangent space and transform data to comply with tg space metric
    sbary <- sqm(bary)
    for (i in 1:size[3]) {
        mydat.tg[, , i] <- sbary %*% Riemanlog(sigma[, , i], bary)
    }

    # Vectorize
    dim(mydat.tg) <- c(size[1] * size[1], size[3])
    # mydat.vec <- matofmat2matofvec(mydat.tg)

    # Performing PCA
    pca.tg <- prcomp(t(mydat.tg))

    # reshape first eigenfunction and pull back
    eigenf <- array(dim = size)
    for (i in 1:size[3]) {
        eigenf[, , i] <- Rieman_exp(matrix(pca.tg$rotation[, i], nrow = size[1], byrow = TRUE), bary)
    }
    pca.tg$eigenf <- eigenf
    pca.tg
}


Riemanlog <- function(A, bary) {
    optGaussMap(bary, A) - diag(dim(A)[1])
}

Rieman_exp <- function(tgB, bary) {
    temp <- tgB + diag(dim(tgB)[1])
    temp %*% bary %*% temp
}
