# Unit tests for NMFN_S4.R
#
# Run from the directory containing NMFN_S4.R:
#   Rscript -e 'library(testthat); source("NMFN_S4.R"); test_file("test_NMFN_S4.R")'
# or inside R:
#   library(testthat); source("NMFN_S4.R"); test_file("test_NMFN_S4.R")
#
# Requires: testthat (install.packages("testthat"))

suppressPackageStartupMessages(library(testthat))
if (!exists("nnmf", mode = "function"))
    suppressMessages(source("NMFN_S4.R"))

# A small, strictly positive, exactly rank-2 matrix: every algorithm should
# be able to reconstruct it well.
set.seed(1)
Wtrue <- matrix(runif(6 * 2), 6, 2)
Htrue <- matrix(runif(2 * 5), 2, 5)
X     <- Wtrue %*% Htrue
Xint  <- matrix(1:12, 3, 4)   # the example from the original docs

rel_err <- function(fit) mean(abs(residuals(fit))) / mean(fit@target)

## ---------------------------------------------------------------------------
context("helpers: mpinv and distance2")

test_that("mpinv inverts a square non-singular matrix", {
    A <- matrix(c(2, 1, 1, 3), 2, 2)
    expect_equal(mpinv(A), solve(A))
})

test_that("mpinv satisfies the Moore-Penrose conditions on a rank-deficient matrix", {
    A  <- matrix(c(1, 2, 2, 4), 2, 2)      # rank 1
    Ap <- mpinv(A)
    expect_equal(A %*% Ap %*% A, A)
    expect_equal(Ap %*% A %*% Ap, Ap)
    expect_equal(t(A %*% Ap), A %*% Ap)
    expect_equal(t(Ap %*% A), Ap %*% A)
})

test_that("mpinv handles a non-square matrix", {
    A  <- matrix(1:6, 2, 3)
    Ap <- mpinv(A)
    expect_equal(dim(Ap), c(3, 2))
    expect_equal(A %*% Ap %*% A, A)
})

test_that("distance2 is the squared Frobenius distance", {
    a <- matrix(1:4, 2); b <- matrix(c(2, 2, 2, 2), 2)
    expect_equal(distance2(a, b), sum((a - b)^2))
    expect_equal(distance2(a, a), 0)
})

## ---------------------------------------------------------------------------
context("algorithm classes")

test_that("constructors return the right classes and labels", {
    expect_s4_class(NMFMultiplicative(), "NMFMultiplicative")
    expect_s4_class(NMFALS(),            "NMFALS")
    expect_s4_class(NMFMultinomial(),    "NMFMultinomial")
    expect_true(is(NMFALS(), "NMFAlgorithm"))
    expect_equal(NMFMultiplicative()@name, "nnmf_mm")
    expect_equal(NMFALS()@name,            "nnmf_als")
    expect_equal(NMFMultinomial()@name,    "nnmf_prob")
})

test_that("NMFAlgorithm is virtual", {
    expect_error(new("NMFAlgorithm"))
})

test_that("nmfAlgorithm maps legacy strings to algorithm objects", {
    expect_s4_class(nmfAlgorithm("nnmf_mm"),   "NMFMultiplicative")
    expect_s4_class(nmfAlgorithm("nnmf_als"),  "NMFALS")
    expect_s4_class(nmfAlgorithm("nnmf_prob"), "NMFMultinomial")
    expect_identical(nmfAlgorithm(NMFALS()), NMFALS())
})

test_that("nmfAlgorithm falls back to multiplicative with a warning", {
    expect_warning(alg <- nmfAlgorithm("bogus"), "unknown method")
    expect_s4_class(alg, "NMFMultiplicative")
})

## ---------------------------------------------------------------------------
context("nnmf(): dispatch and return type")

test_that("default method is multiplicative updates", {
    set.seed(10)
    fit <- nnmf(Xint, 2, verbose = FALSE)
    expect_s4_class(fit, "NMFfit")
    expect_s4_class(algorithm(fit), "NMFMultiplicative")
})

test_that("method strings dispatch to the right algorithm", {
    set.seed(11)
    expect_s4_class(algorithm(nnmf(Xint, 2, "nnmf_mm",   verbose = FALSE)), "NMFMultiplicative")
    expect_s4_class(algorithm(nnmf(Xint, 2, "nnmf_als",  verbose = FALSE)), "NMFALS")
    expect_s4_class(algorithm(nnmf(Xint, 2, "nnmf_prob", verbose = FALSE)), "NMFMultinomial")
})

test_that("algorithm objects are accepted as 'method'", {
    set.seed(12)
    fit <- nnmf(Xint, 2, NMFALS(), verbose = FALSE)
    expect_s4_class(algorithm(fit), "NMFALS")
})

test_that("nmfFit() can be called directly", {
    set.seed(13)
    fit <- nmfFit(NMFMultinomial(), Xint, 2, maxiter = 100, verbose = FALSE)
    expect_s4_class(fit, "NMFfit")
    expect_s4_class(algorithm(fit), "NMFMultinomial")
})

test_that("verbose = TRUE prints the algorithm label and progress", {
    set.seed(14)
    out <- capture.output(fit <- nnmf(Xint, 2, "nnmf_als", maxiter = 100))
    expect_match(out[1], "Alternating Least Squares Algorithm")
    expect_true(any(grepl("Iter = ", out)))
})

test_that("verbose = FALSE prints nothing", {
    set.seed(15)
    out <- capture.output(fit <- nnmf(Xint, 2, maxiter = 100, verbose = FALSE))
    expect_length(out, 0)
})

# run a legacy wrapper with its printed progress swallowed
capture_fit <- function(f, ...) { fit <- NULL; capture.output(fit <- f(...)); fit }

test_that("legacy wrapper functions still exist and return NMFfit", {
    set.seed(16)
    expect_s4_class(capture_fit(nnmf_mm,   Xint, 2, 100), "NMFfit")
    expect_s4_class(capture_fit(nnmf_als,  Xint, 2, 100), "NMFfit")
    expect_s4_class(capture_fit(nnmf_prob, Xint, 2, 100), "NMFfit")
})

## ---------------------------------------------------------------------------
context("NMFfit: shape, validity and accessors")

for (alg in list(NMFMultiplicative(), NMFALS(), NMFMultinomial())) {
    test_that(paste0(alg@name, ": W and H have the right shape and are non-negative"), {
        set.seed(20)
        fit <- nnmf(X, 2, alg, maxiter = 200, verbose = FALSE)
        expect_equal(dim(basis(fit)), c(6, 2))
        expect_equal(dim(coef(fit)),  c(2, 5))
        expect_true(all(basis(fit) >= 0))
        expect_true(all(coef(fit)  >= 0))
        expect_true(validObject(fit))
        expect_equal(dim(fit), dim(X))
        expect_equal(nfactors(fit), 2)
    })
}

test_that("fitted() and residuals() are consistent with W and H", {
    set.seed(21)
    fit <- nnmf(X, 2, maxiter = 100, verbose = FALSE)
    expect_equal(fitted(fit), basis(fit) %*% coef(fit))
    expect_equal(residuals(fit), X - fitted(fit))
    expect_equal(fit@target, X)
})

test_that("legacy $W / $H and names() work", {
    set.seed(22)
    fit <- nnmf(Xint, 2, maxiter = 100, verbose = FALSE)
    expect_identical(fit$W, basis(fit))
    expect_identical(fit$H, coef(fit))
    expect_equal(names(fit), c("W", "H"))
    expect_error(fit$foo, "no component")
})

test_that("niter() and converged() are populated", {
    set.seed(23)
    fit <- nnmf(X, 2, maxiter = 100, verbose = FALSE)
    expect_type(niter(fit), "integer")
    expect_true(niter(fit) >= 1 && niter(fit) <= 100)
    expect_type(converged(fit), "logical")
    expect_length(converged(fit), 1)
})

test_that("the trace records one row per 50 iterations", {
    set.seed(24)
    fit <- nnmf(X, 2, maxiter = 150, verbose = FALSE)
    tr <- nmfTrace(fit)
    expect_s3_class(tr, "data.frame")
    expect_named(tr, c("iteration", "relError", "diff", "euclDist"))
    expect_equal(nrow(tr), niter(fit) %/% 50)
    expect_equal(tr$iteration, seq(50, niter(fit), by = 50))
})

test_that("show() and summary() run and mention the algorithm", {
    set.seed(25)
    fit <- nnmf(Xint, 2, "nnmf_als", maxiter = 100, verbose = FALSE)
    out <- capture.output(show(fit))
    expect_true(any(grepl("NMFfit", out)))
    expect_true(any(grepl("Alternating Least Squares", out)))
    out <- capture.output(s <- summary(fit))
    expect_named(s, c("W", "H", "trace"))
})

test_that("the call is stored", {
    set.seed(26)
    fit <- nnmf(Xint, 2, maxiter = 100, verbose = FALSE)
    expect_type(fit@call, "language")
    expect_equal(as.character(fit@call[[1]]), "nnmf")
})

## ---------------------------------------------------------------------------
context("validity")

test_that("validity rejects mismatched dimensions and negative entries", {
    set.seed(27)
    fit <- nnmf(X, 2, maxiter = 50, verbose = FALSE)
    bad <- fit; bad@H <- bad@H[1, , drop = FALSE]
    expect_error(validObject(bad), "ncol\\(W\\)")
    bad <- fit; bad@W[1, 1] <- -1
    expect_error(validObject(bad), "non-negative")
    bad <- fit; bad@target <- bad@target[1:3, ]
    expect_error(validObject(bad), "dim\\(W")
})

## ---------------------------------------------------------------------------
context("numerical behaviour")

for (alg in list(NMFMultiplicative(), NMFALS(), NMFMultinomial())) {
    test_that(paste0(alg@name, ": reconstructs an exactly rank-2 matrix"), {
        set.seed(30)
        fit <- nnmf(X, 2, alg, maxiter = 1000, verbose = FALSE)
        expect_lt(rel_err(fit), 0.05)
    })
}

test_that("more factors do not fit worse (multiplicative)", {
    set.seed(31)
    e1 <- rel_err(nnmf(X, 1, maxiter = 500, verbose = FALSE))
    set.seed(31)
    e3 <- rel_err(nnmf(X, 3, maxiter = 500, verbose = FALSE))
    expect_lte(e3, e1 + 1e-8)
})

test_that("results are reproducible under set.seed", {
    set.seed(32); a <- nnmf(X, 2, "nnmf_als", maxiter = 100, verbose = FALSE)
    set.seed(32); b <- nnmf(X, 2, "nnmf_als", maxiter = 100, verbose = FALSE)
    expect_equal(basis(a), basis(b))
    expect_equal(coef(a),  coef(b))
})

test_that("maxiter is respected", {
    set.seed(33)
    fit <- nnmf(X, 2, maxiter = 7, verbose = FALSE)
    expect_equal(niter(fit), 7L)
    expect_false(converged(fit))
    expect_equal(nrow(nmfTrace(fit)), 0)
})

test_that("data.frame input is accepted", {
    set.seed(34)
    fit <- nnmf(as.data.frame(X), 2, maxiter = 50, verbose = FALSE)
    expect_s4_class(fit, "NMFfit")
    expect_equal(dim(fit), dim(X))
})

## ---------------------------------------------------------------------------
context("input validation")

test_that("non-finite input is rejected", {
    bad <- X; bad[1, 1] <- NA
    expect_error(nnmf(bad, 2, verbose = FALSE), "infinite or missing")
    bad[1, 1] <- Inf
    expect_error(nnmf(bad, 2, verbose = FALSE), "infinite or missing")
})

test_that("empty and non-numeric input is rejected", {
    expect_error(nnmf(matrix(numeric(0), 0, 3), 2, verbose = FALSE), "0 extent")
    expect_error(nnmf(matrix(letters[1:4], 2), 2, verbose = FALSE), "numeric")
})

test_that("invalid k is rejected", {
    expect_error(nnmf(X, 0,    verbose = FALSE), "'k'")
    expect_error(nnmf(X, 1.5,  verbose = FALSE), "'k'")
    expect_error(nnmf(X, -2,   verbose = FALSE), "'k'")
    expect_error(nnmf(X, c(1, 2), verbose = FALSE), "'k'")
})
