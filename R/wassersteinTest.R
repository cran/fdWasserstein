wassersteinTest <- function(data, grp, B = 1000,
                            statistic = c("transport", "distance"),
                            type = c("permutation", "bootstrap"),
                            r = c("HS", "trace", "operator"),
                            align = TRUE,
                            use.future = FALSE,
                            iter.bary = 10) {
    statistic <- match.arg(statistic)
    r <- match.arg(r)
    rep <- (type <- match.arg(type)) == "bootstrap"
    level <- unique(grp)
    g <- sapply(level, function(l) grp == l)
    freq <- colSums(g)

    if (align) {
        for (i in seq_along(level)) data[g[, i], ] <- scale(data[g[, i], ], scale = FALSE)
    }
    d <- c(NCOL(data), NCOL(data), length(level))
    Test <- if (statistic == "transport") {
        Id <- diag(d[1])
        ff <- rep(freq, rep(d[1] * d[1], d[3]))
        function(idx) {
            S <- array(apply(g, 2, function(gi) cov(data[idx[gi], ])), d)
            Sbar <- gaussBary(S, freq, max.iter = iter.bary, silent = TRUE)$gamma
            T <- apply(S, 3, function(Si) optGaussMap(Sbar, Si) - Id)
            if (r == "HS"){
                sum(ff * T^2)
            }else if(r == "trace"){
                #sum(abs(.svd(A)))^2
                sum(ff*(sum(diag(.svd(T))))^2)
            (nuclear.norm(T))^(2)
            }else if (r == "operator"){
            sum(ff*(operator.norm(T))^(2))
            }
        }
    } else {
        function(idx) {
            S <- array(apply(g, 2, function(gi) cov(data[idx[gi], ])), d)
            Sbar <- gaussBary(S, freq, max.iter = iter.bary, silent = TRUE)$gamma
            sum(freq * apply(S, 3, dwstein, Sbar))
        }
    }
    idx <- 1:NROW(data)
    obs <- Test(idx)
    if (use.future) {
        u <- list()
        for (i in seq.int(B)) u[[i]] <- future::future(Test(sample(idx, replace = rep)), seed = TRUE)
        trep <- sapply(u, future::value)
    } else {
        trep <- replicate(B, Test(sample(idx, replace = rep)))
    }
    list(stat = obs, p.value = (1 + sum(trep > obs)) / (1 + B), trep = trep)
}




.svd <-function(A){
    ATA <- t(A) %*% A
    return(sqm(ATA))
    #return(sqrt(.eigen(ATA)$values))
}

nuclear.norm <- function(A){
    sum(abs(.svd(A)))^2
}


operator.norm <- function (x) 
{
    if (!is.numeric(x)) {
        stop("argument x is not numeric")
    }
    if (is.vector(x)) {
        return(sqm(sum(x * x)))
    }
    if (!is.matrix(x)) {
        return("argument x is not a matrix")
    }
    A <- t(x) %*% x
    eigenA <- .eigen(A)
    lambdaA <- eigenA$values
    maxLambdaA <- lambdaA[1]
    if (maxLambdaA < 0) {
        stop("t(x) %*% x is negative definite")
    }
    return(maxLambdaA)
}
