# Non-negative Matrix Factorization (NMFN) -- S4 rewrite
#
# Algorithms: multiplicative updates (Lee & Seung 2001), alternating least
# squares (after Lars Kai Hansen's nnmf_als Matlab code) and a multinomial
# algorithm. Moore-Penrose inverse after Torsten Hothorn.
# Original R code by Timothy Liu; this is a rewrite of that code in S4.
#
# Design
#   * Each algorithm is an S4 class deriving from the virtual class
#     "NMFAlgorithm". The fitting code is a method of the generic nmfFit(),
#     so the algorithm is chosen by S4 dispatch instead of string comparison.
#     Adding an algorithm = one setClass() + one setMethod("nmfFit", ...).
#   * A fit is returned as an "NMFfit" object holding W, H, the target matrix,
#     the algorithm used, the iteration count, a convergence flag and a trace
#     of the progress that the original code printed with cat().
#   * nnmf() keeps the original signature and still accepts the method strings
#     "nnmf_mm", "nnmf_als", "nnmf_prob"; it also accepts an algorithm object.
#   * fit$W / fit$H still work, so old code that did z$W keeps running.
#
# Usage
#   source("NMFN_S4.R")
#   X   <- matrix(1:12, 3, 4)
#   fit <- nnmf(X, 2)                     # multiplicative update (default)
#   fit <- nnmf(X, 2, "nnmf_als")         # alternating least squares
#   fit <- nnmf(X, 2, NMFMultinomial())   # multinomial, via algorithm object
#   fit                                   # show()
#   basis(fit); coef(fit); fitted(fit); residuals(fit); summary(fit)

suppressPackageStartupMessages(library(methods))

## ===========================================================================
## 1. Helper functions (unchanged numerics from the original)
## ===========================================================================

# Moore-Penrose inverse
mpinv <- function(X)
{
    Eps <- 100 * .Machine$double.eps
    s <- svd(X)
    d <- s$d
    m <- length(d)

    if (!(is.vector(d)))
        return(t(s$v %*% (1/d) %*% t(s$u)))

    d <- d[d > Eps]
    notnull <- length(d)

    if (notnull == 1) {
        inv <- 1/d
    } else {
        inv <- solve(diag(d))
    }

    if (notnull != m) {
        inv <- cbind(inv, matrix(0, nrow = notnull, ncol = (m - notnull)))
        inv <- rbind(inv, matrix(0, nrow = (m - notnull), ncol = m))
    }

    mp <- s$v %*% inv %*% t(s$u)
    mp[abs(mp) < Eps] <- 0
    mp
}

# Squared Euclidean distance between two matrices
distance2 <- function(x1, x2)
{
    temp <- x1 - x2
    sum(temp * temp)
}

.PRINT_ITER <- 50L    # iterations between progress checks (as in the original)
.TOL        <- 10e-6  # relative-error stopping tolerance (as in the original)

.checkInput <- function(x, k)
{
    x <- as.matrix(x)
    if (!is.numeric(x))
        stop("'x' must be numeric")
    if (any(!is.finite(x)))
        stop("infinite or missing values in 'x'")
    if (!nrow(x) || !ncol(x))
        stop("0 extent dimensions")
    if (length(k) != 1L || !is.numeric(k) || k < 1 || k != round(k))
        stop("'k' must be a single positive integer")
    storage.mode(x) <- "double"
    x
}

.emptyTrace <- function()
    data.frame(iteration = integer(), relError = numeric(),
               diff = numeric(), euclDist = numeric())

# One progress check (what the original did inside `if (iter %% 50 == 0)`).
# Returns the updated reconstruction, the trace and whether to stop.
.progress <- function(iter, x, W, H, Xr_old, trace, verbose)
{
    Xr   <- W %*% H
    diff <- sum(abs(Xr_old - Xr))
    eucl <- distance2(x, Xr)
    err  <- mean(abs(x - Xr)) / mean(x)
    if (verbose)
        cat("Iter = ", iter, "\t", "relative error = ", err, "\t",
            "diff = ", diff, "\t", "eucl dist = ", eucl, "\n")
    trace <- rbind(trace, data.frame(iteration = iter, relError = err,
                                     diff = diff, euclDist = eucl))
    stop_now <- err < .TOL
    if (stop_now && verbose)
        cat("Execution finishes at iteration = ", iter, "\n")
    list(Xr_old = Xr, trace = trace, stop = stop_now)
}

## ===========================================================================
## 2. Classes
## ===========================================================================

# --- algorithm classes ------------------------------------------------------
setClass("NMFAlgorithm",
         representation("VIRTUAL",
                        name  = "character",   # legacy method string
                        label = "character"))  # human-readable name

setClass("NMFMultiplicative", contains = "NMFAlgorithm",
         prototype = prototype(name = "nnmf_mm",
                   label = "Multiplicative Update Algorithm"))

setClass("NMFALS", contains = "NMFAlgorithm",
         prototype = prototype(name = "nnmf_als",
                   label = "Alternating Least Squares Algorithm"))

setClass("NMFMultinomial", contains = "NMFAlgorithm",
         prototype = prototype(name = "nnmf_prob",
                   label = "Multinomial Algorithm"))

# Constructors
NMFMultiplicative <- function() new("NMFMultiplicative")
NMFALS            <- function() new("NMFALS")
NMFMultinomial    <- function() new("NMFMultinomial")

# Map a legacy method string to an algorithm object.
# Unknown strings fall back to multiplicative updates, as the original did.
nmfAlgorithm <- function(method = "nnmf_mm")
{
    if (is(method, "NMFAlgorithm")) return(method)
    switch(as.character(method),
           nnmf_als  = NMFALS(),
           nnmf_prob = NMFMultinomial(),
           nnmf_mm   = NMFMultiplicative(),
           {
               warning("unknown method '", method,
                       "'; using multiplicative updates")
               NMFMultiplicative()
           })
}

# --- result class -----------------------------------------------------------
setClass("NMFfit",
         representation(W          = "matrix",      # n x k basis
                        H          = "matrix",      # k x m coefficients
                        algorithm  = "NMFAlgorithm",
                        target     = "matrix",      # matrix that was factorised
                        iterations = "integer",
                        converged  = "logical",
                        trace      = "data.frame",  # progress every 50 iters
                        call       = "call"))

setValidity("NMFfit", function(object) {
    msgs <- character()
    W <- object@W; H <- object@H; X <- object@target
    if (ncol(W) != nrow(H))
        msgs <- c(msgs, sprintf("ncol(W) (%d) must equal nrow(H) (%d)",
                                ncol(W), nrow(H)))
    if (length(X) && (nrow(W) != nrow(X) || ncol(H) != ncol(X)))
        msgs <- c(msgs, "dim(W %*% H) must equal dim(target)")
    if (any(W < 0, na.rm = TRUE) || any(H < 0, na.rm = TRUE))
        msgs <- c(msgs, "W and H must be non-negative")
    if (length(msgs)) msgs else TRUE
})

# Internal constructor used by the fitting methods
.NMFfit <- function(W, H, algorithm, x, iterations, converged, trace, call)
{
    # clamp tiny negative round-off (e.g. -1e-17) so validity holds
    W[W < 0 & W > -1e-12] <- 0
    H[H < 0 & H > -1e-12] <- 0
    new("NMFfit", W = W, H = H, algorithm = algorithm, target = x,
        iterations = as.integer(iterations), converged = converged,
        trace = trace, call = call)
}

## ===========================================================================
## 3. Generics
## ===========================================================================

setGeneric("nnmf",
           function(x, k, method = "nnmf_mm", maxiter = 1000,
                    eps = 2.2204e-16, verbose = TRUE, ...)
               standardGeneric("nnmf"),
           signature = c("x", "method"))

setGeneric("nmfFit",
           function(algorithm, x, k, maxiter = 1000, eps = 2.2204e-16,
                    verbose = TRUE, ...)
               standardGeneric("nmfFit"),
           signature = "algorithm")

setGeneric("basis",     function(object, ...) standardGeneric("basis"))
setGeneric("algorithm", function(object, ...) standardGeneric("algorithm"))
setGeneric("nfactors",  function(object, ...) standardGeneric("nfactors"))
setGeneric("niter",     function(object, ...) standardGeneric("niter"))
setGeneric("converged", function(object, ...) standardGeneric("converged"))
setGeneric("nmfTrace",  function(object, ...) standardGeneric("nmfTrace"))

## ===========================================================================
## 4. Fitting methods -- one per algorithm class
## ===========================================================================

# --- multiplicative updates (Lee & Seung 2001) ------------------------------
setMethod("nmfFit", "NMFMultiplicative",
function(algorithm, x, k, maxiter = 1000, eps = 2.2204e-16, verbose = TRUE, ...)
{
    cl <- match.call()
    x  <- .checkInput(x, k)
    n  <- nrow(x); m <- ncol(x)

    W <- matrix(abs(rnorm(n * k)), n, k)
    H <- matrix(abs(rnorm(k * m)), k, m)

    Xr_old <- W %*% H
    trace  <- .emptyTrace()
    converged <- FALSE
    iter <- 0L

    for (iter in seq_len(maxiter)) {
        H <- H * (t(W) %*% x) / ((t(W) %*% W) %*% H + eps)
        W <- W * t(H %*% t(x)) / (W %*% (H %*% t(H)) + eps)

        if (iter %% .PRINT_ITER == 0) {
            p <- .progress(iter, x, W, H, Xr_old, trace, verbose)
            Xr_old <- p$Xr_old; trace <- p$trace
            if (p$stop) { converged <- TRUE; break }
        }
    }
    .NMFfit(W, H, algorithm, x, iter, converged, trace, cl)
})

# --- alternating least squares ---------------------------------------------
setMethod("nmfFit", "NMFALS",
function(algorithm, x, k, maxiter = 1000, eps = 2.2204e-16, verbose = TRUE, ...)
{
    cl <- match.call()
    x  <- .checkInput(x, k)
    D  <- nrow(x); N <- ncol(x)
    Xscale <- sum(x)

    W <- matrix(abs(rnorm(D * k)), D, k)
    H <- matrix(abs(rnorm(k * N)), k, N)

    Rscale  <- sum(W %*% H)
    sqrnorm <- sqrt(Rscale / Xscale)
    H <- H / sqrnorm
    W <- W / sqrnorm

    Xr_old <- W %*% H
    trace  <- .emptyTrace()
    converged <- FALSE
    iter <- 0L

    for (iter in seq_len(maxiter)) {
        W <- x %*% t(mpinv(H %*% t(H)) %*% H)
        W <- (W > 0) * W
        W <- W / (t(matrix(rep(colSums(W), D), ncol(W), nrow(W))) + eps)

        H <- t(W %*% mpinv(t(W) %*% W)) %*% x
        H <- H * (H > 0)

        if (iter %% .PRINT_ITER == 0) {
            p <- .progress(iter, x, W, H, Xr_old, trace, verbose)
            Xr_old <- p$Xr_old; trace <- p$trace
            if (p$stop) { converged <- TRUE; break }
        }
    }
    .NMFfit(W, H, algorithm, x, iter, converged, trace, cl)
})

# --- multinomial -------------------------------------------------------------
setMethod("nmfFit", "NMFMultinomial",
function(algorithm, x, k, maxiter = 1000, eps = 100 * .Machine$double.eps,
         verbose = TRUE, ...)
{
    cl <- match.call()
    x  <- .checkInput(x, k)
    D  <- nrow(x); N <- ncol(x)

    X_factor <- sum(x)
    X_org    <- x
    x <- x / X_factor

    W <- matrix(abs(rnorm(D * k)), D, k)
    W <- W / t(matrix(rep(colSums(W), D), ncol(W), nrow(W)))

    H <- matrix(abs(rnorm(k * N)), k, N)
    H <- H / (matrix(rep(rowSums(H), N), nrow(H), ncol(H)))

    P <- matrix(1, k, 1)
    P <- P / sum(P)

    W1 <- W; H1 <- H
    Xr_old <- W %*% H
    trace  <- .emptyTrace()
    converged <- FALSE
    iter <- 0L

    for (iter in seq_len(maxiter)) {
        Qnorm <- (W %*% diag(diag(P), k, k)) %*% H

        for (j in seq_len(k)) {
            Q  <- (t(t(W[, j])) %*% H[j, ] * P[j]) / (Qnorm + eps)
            XQ <- x * Q

            dummy   <- rowSums(XQ)
            W1[, j] <- t(dummy / sum(dummy))

            dummy   <- colSums(XQ)
            H1[j, ] <- (dummy / sum(dummy))
        }
        W <- W1; H <- H1

        if (iter %% .PRINT_ITER == 0) {
            p <- .progress(iter, x, W, H, Xr_old, trace, verbose)
            Xr_old <- p$Xr_old; trace <- p$trace
            if (p$stop) { converged <- TRUE; break }
        }
    }

    W <- W %*% diag(diag(sqrt(P)), k, k) * X_factor
    H <- diag(diag(sqrt(P)), k, k) %*% H
    .NMFfit(W, H, algorithm, X_org, iter, converged, trace, cl)
})

## ===========================================================================
## 5. nnmf() -- the user-facing entry point (same signature as the original)
## ===========================================================================

# method given as a string ("nnmf_mm", "nnmf_als", "nnmf_prob")
setMethod("nnmf", signature(x = "ANY", method = "character"),
function(x, k, method = "nnmf_mm", maxiter = 1000, eps = 2.2204e-16,
         verbose = TRUE, ...)
{
    alg <- nmfAlgorithm(method)
    callGeneric(x = x, k = k, method = alg, maxiter = maxiter, eps = eps,
                verbose = verbose, ...)
})

# method given as an algorithm object
setMethod("nnmf", signature(x = "ANY", method = "NMFAlgorithm"),
function(x, k, method = "nnmf_mm", maxiter = 1000, eps = 2.2204e-16,
         verbose = TRUE, ...)
{
    if (verbose) cat(method@label, "\n")
    fit <- nmfFit(method, x, k, maxiter = maxiter, eps = eps,
                  verbose = verbose, ...)
    fit@call <- match.call()
    fit
})

# method missing -> default algorithm
setMethod("nnmf", signature(x = "ANY", method = "missing"),
function(x, k, method, maxiter = 1000, eps = 2.2204e-16, verbose = TRUE, ...)
    nnmf(x, k, NMFMultiplicative(), maxiter = maxiter, eps = eps,
         verbose = verbose, ...))

# Backward-compatible wrappers with the original names and signatures
nnmf_mm   <- function(x, k, maxiter = 1000, eps = 2.2204e-16)
    nmfFit(NMFMultiplicative(), x, k, maxiter, eps)
nnmf_als  <- function(x, k, maxiter = 1000, eps = 2.2204e-16)
    nmfFit(NMFALS(), x, k, maxiter, eps)
nnmf_prob <- function(x, k, maxiter = 1000, eps = 100 * .Machine$double.eps)
    nmfFit(NMFMultinomial(), x, k, maxiter, eps)

## ===========================================================================
## 6. Methods on NMFfit
## ===========================================================================

setMethod("basis",     "NMFfit", function(object, ...) object@W)
setMethod("coef",      "NMFfit", function(object, ...) object@H)
setMethod("algorithm", "NMFfit", function(object, ...) object@algorithm)
setMethod("nfactors",  "NMFfit", function(object, ...) ncol(object@W))
setMethod("niter",     "NMFfit", function(object, ...) object@iterations)
setMethod("converged", "NMFfit", function(object, ...) object@converged)
setMethod("nmfTrace",  "NMFfit", function(object, ...) object@trace)
setMethod("fitted",    "NMFfit", function(object, ...) object@W %*% object@H)
setMethod("residuals", "NMFfit",
          function(object, ...) object@target - object@W %*% object@H)
setMethod("dim",       "NMFfit", function(x) dim(x@target))

# Legacy list-style access: fit$W, fit$H, names(fit)
setMethod("names", "NMFfit", function(x) c("W", "H"))
setMethod("$", "NMFfit", function(x, name)
{
    switch(name,
           W = x@W,
           H = x@H,
           stop("no component '", name, "'; use W or H, or the accessors ",
                "basis()/coef()"))
})

setMethod("show", "NMFfit", function(object)
{
    d <- dim(object)
    cat("<Object of class NMFfit>\n")
    cat("  algorithm : ", object@algorithm@label,
        " (", object@algorithm@name, ")\n", sep = "")
    cat("  target    : ", d[1], " x ", d[2], " matrix\n", sep = "")
    cat("  factors   : ", nfactors(object), "\n", sep = "")
    cat("  W         : ", nrow(object@W), " x ", ncol(object@W), "\n", sep = "")
    cat("  H         : ", nrow(object@H), " x ", ncol(object@H), "\n", sep = "")
    cat("  iterations: ", object@iterations,
        if (object@converged) " (converged)" else " (maxiter reached)",
        "\n", sep = "")
    err <- mean(abs(residuals(object))) / mean(object@target)
    cat("  rel. error: ", format(err, digits = 4), "\n", sep = "")
    invisible(object)
})

setMethod("summary", "NMFfit", function(object, ...)
{
    show(object)
    tr <- object@trace
    if (nrow(tr)) {
        cat("\nProgress (every", .PRINT_ITER, "iterations):\n")
        print(tr, row.names = FALSE)
    }
    invisible(list(W = object@W, H = object@H, trace = tr))
})
