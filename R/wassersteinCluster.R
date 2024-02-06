wassersteinCluster <- function(data, grp, kmin = 2, kmax = 10,
                               E = -0.75 * (0.95 * log(0.95) + 0.05 * log(0.05)) + 0.25 * log(2),
                               nstart = 5, nrefine = 5, ntry=0,
                               max.iter = 20, tol = 0.001,
                               nreduced = length(unique(grp)), 
                               nperm=0,
                               add.sigma=FALSE, 
                               use.future=FALSE, 
                               verbose = TRUE) {
    if (use.future && !requireNamespace("future", quietly=TRUE)) 
        stop("please, install package future")
    if (min(table(grp))<2) stop("a minimum of two curves for each sample is needed")
    a <- .wassersteinCluster(data, grp, kmin, kmax, E, nstart, nrefine, ntry, max.iter, tol, nreduced, 
                             add.sigma, use.future, verbose)
    if (nperm>0) {
        if (verbose) cat("Permutation test - please wait!\n")
        for (i in unique(grp)) data[grp==i,] <- scale(data[grp==i,], scale=FALSE)
        fone <- function(g) 
            max(trimmedAverageSilhouette(
                .wassersteinCluster(data, g, kmin, kmax, 
                                    E, nstart, nrefine, ntry, max.iter, tol, 
                                    nreduced,FALSE, FALSE, FALSE),
                FALSE
          ))
        if (use.future) {
            if (verbose) cat("Using package future...")
            u <- list()
            for (i in seq.int(nperm)) u[[i]] <- future::future(fone(sample(grp)), seed=TRUE)
            tasw0 <-  sapply(u, future::value)
        } else {
            tasw0 <- numeric(nperm)
            for (i in seq.int(nperm)) {
                if (verbose) cat(".")
                tasw0[i] <- fone(sample(grp))
            }
        }
        if (verbose) cat("\n")
        v <- max(trimmedAverageSilhouette(a,FALSE))
        attr(a, "tasw.test") <- list(
                                     maxTASW=v,
                                     tasw0=tasw0, 
                                     p.value=(1+sum(v<tasw0))/(nperm+2)
                                     )
    }
    a
}

.wassersteinCluster <- function(X, g, kmin, kmax, E, nstart, nrefine, ntry, max.iter, tol, nreduced, 
                                add.sigma, use.future, verbose) {
    g <- as.factor(g)
    lg <- levels(g)
    N <- length(lg)
    M <- NCOL(X)
    sigma <- array(NA, c(M, M, N))
    df <- numeric(N)
    for (i in seq_along(lg)) {
        jdx <- which(g == lg[i])
        sigma[,,i] <- cov(X[jdx,])
        df[i] <- length(jdx) - 1
    }
    if (nreduced<N) {
        idx <- sample(seq.int(N), nreduced)
        fsigma <- sigma
        fdf <- df
        sigma <- sigma[,,idx]
        df <- df[idx]
    }
    sigma0.5 <- array(apply(sigma, 3, sqm), dim(sigma))
    dw <- local({
        cache <- new.env(hash = TRUE)
        function(i, j) {
            if (i > j) {
                tmp <- i
                i <- j
                j <- tmp
            }
            nm <- as.character(i * N + j)
            if (is.null(r <- get0(nm, cache))) {
                assign(nm, dwstein(sigma[, , i], sigma[, , j], sigma0.5[, , j]), cache)
            } else {
                r
            }
        }
    })
    ans <- vector("list", kmax - kmin + 1)
    lone <- function(init) {
        wassersteinClusterLocal(sigma, df, sigma0.5, sigma[,,init$ibest], init$w, init$obest, 
                                E, max.iter, tol, verbose)
    }
    for (k in kmin:kmax) {
       init <- wassersteinClusterInit(sigma, df, k, E, nstart, nrefine, ntry, verbose, dw)
       ans[[k-kmin+1]] <- if (use.future) future::future(lone(init)) else lone(init)
    }
    if (use.future) ans <- lapply(ans, future::value)
    if (nreduced < N) {
        ans <- lapply(ans,  
                      function(a){
                        a$d <- apply(a$g, 3, function(g) apply(fsigma, 3, dwstein, g, sqm(g)))
                        a$w <- findWeights(fdf*a$d, a$E)
                        a$obj <- sum(fdf*a$w*a$d)
                        a
                    })
    }
    if (add.sigma) {
        attr(ans, "sample.covariances") <- if (nreduced < N) fsigma else sigma
    }
    attr(ans, "df") <- if (nreduced < N) fdf else df
    ans
}

trimmedAverageSilhouette <- function(a, plot=TRUE) {
    df <- attr(a, "df")
    tasw <- sapply(a, 
                   function(b) {
                      C <- apply(b$w, 1, max)
                      ic <- which(C>=mean(C))
                      dd <- apply(b$d, 1, sort)[1:2,]
                      s <- 1-dd[1,]/dd[2,]
                      sum(s[ic]*df[ic]) / sum(df[ic])
                   })
    nc <- sapply(a, function(b) b$k)
    names(tasw) <- nc
    if (plot) {
       plot(nc, tasw, type="b", 
            ylab="trimmed average silhouette width", 
            xlab="number of clusters",
            xaxt="n")
            axis(1, at=nc, labels=nc)
    }
    tasw
}


wassersteinClusterInit <- function(sigma, df, k, E, nstart, nrefine, ntry, verbose, dw) {
    obest <- +Inf
    n <- dim(sigma)[3]
    sn <- seq_len(n)
    d <- matrix(0, n, k)
    idx <- integer(k)
    sk <- seq_len(k)
    ntry <- if (ntry==0) round(1 + n / k) else ntry
    if (verbose) cat("Initialization k:", k, "\n")
    for (ns in seq_len(nstart)) {
        p <- rep(1, n)
        for (i in sk) {
            idx[i] <- j <- sample(sn, 1, prob = p)
            d[, i] <- sapply(sn, dw, j)
            p <- if (i == 1) d[, i] else pmin(p, d[, i])
        }
        w <- findWeights(df * d, E)
        obj <- sum(df * w * d)
        for (nr in seq_len(nrefine)) {
            for (i in sk) {
                p <- apply(d[, -i, drop = FALSE], 1, min)
                for (j in sample(sn, min(ntry, sum(p > 0)), prob = p)) {
                    dold <- d[, i]
                    d[, i] <- sapply(sn, dw, j)
                    neww <- findWeights(df*d, E)
                    newobj <- sum(df * neww * d)
                    if (newobj < obj) {
                        idx[i] <- j
                        obj <- newobj
                        w <- neww
                    } else {
                        d[, i] <- dold
                    }
                }
            }
        }
        if (obj < obest) {
            obest <- obj
            ibest <- idx
            dbest <- d
            wbest <- w
        }
        if (verbose) {
            cat(
                "    Round:", ns, "Best criterion:", obest,
                " Best subset:", sort(ibest), "\n"
            )
        }
    }
    list(ibest=ibest, w=wbest, obest=obest)
}

wassersteinClusterLocal <- function(sigma, df, sigma0.5, g, w, obest, E, max.iter, tol, verbose) {
    ## Main optimization loop
    k <- NCOL(w)
    n <- NROW(w)
    if (verbose) cat("Optimization k:", k, "\n")
    for (ns in seq_len(max.iter)) {
        for (i in seq_len(k)) {
            g[, , i] <- gaussBary(sigma, df*w[, i], g[,,i], sigma0.5,  max.iter=1, silent=TRUE)$gamma 
        }
        d <- apply(g, 3, function(x) apply(sigma, 3, dwstein, x, sqm(x)))
        w <- findWeights(df*d, E)
        old <- obest
        obest <- sum(df * w * d)
        if (verbose) cat("    Iteration:", ns, "Criterion:", obest, "\n")
        if (old - obest < tol * old) break
    }
    list(
        k = k, E = E, eta = attr(w, "eta"),
        w = w, g = g, d = d,
        obj = obest
    )
}


findWeights <- function(d, E) {
    E <- NROW(d) * E
    d <- d - apply(d, 1, min)
    f <- function(eta) {
        w <- exp(-d / eta)
        w <- w / rowSums(w)
        a <- sum(w * log(w)) + E
        if (is.na(a)) E else a
    }
    eta <- uniroot(f, c(0, 100 * max(d)))$root
    w <- exp(-d / eta)
    w <- w / rowSums(w)
    attr(w, "eta") <- eta
    w
}


