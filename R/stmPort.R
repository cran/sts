##############################################################################
## Internal implementation of the Structural Topic Model (STM) routines that
## sts uses for initialization.
##
## sts previously called stm::stm() to (i) compute the spectral (anchor word)
## initialization of the topic-word distributions (initialization = "anchor")
## and (ii) fit a full STM with a content covariate (initialization = "stm").
## The functions below reproduce exactly the computations performed by
## stm::stm() for those two calls, so that sts no longer depends on the stm
## package. They are adapted from stm version 1.3.8 (files R/stm.R,
## R/stm.control.R, R/STMinit.R, R/spectral.R, R/STMestep.R, R/STMoptbeta.R,
## R/STMmnreg.R, R/STMconvergence.R and R/STMfunctions.R), written by
## Margaret E. Roberts, Brandon M. Stewart and Dustin Tingley and released
## under the MIT license:
##   Copyright (c) 2017 Margaret Roberts, Brandon Stewart and Dustin Tingley
## See inst/COPYRIGHTS for the full notice.
##############################################################################

## Fit (or initialize) an STM model.
##
## Equivalent to
##   stm::stm(documents, vocab, K, prevalence = ~X, data = data.frame(X),
##            content = content, max.em.its = max.em.its, LDAbeta = FALSE,
##            verbose = FALSE)
## with all other arguments at their stm defaults (spectral initialization,
## pooled prevalence prior, L1 content-covariate prior, interactions = TRUE).
## Only the elements of the fitted object that sts uses are returned, using
## the same names and structure as an STM object.
stmInit <- function(documents, vocab, K, X, content = NULL, max.em.its = 0,
                    emtol = 1e-5) {
  ## ---- checks on the documents (as in stm::stm) ----
  if (!is.list(documents)) stop("documents must be a list, see documentation.")
  if (!all(unlist(lapply(documents, is.matrix)))) {
    stop("Each list element in documents must be a matrix. See documentation.")
  }
  if (any(unlist(lapply(documents, function(x) anyDuplicated(x[1, ]))))) {
    stop("Duplicate term indices within a document.  See documentation for proper format.")
  }
  N <- length(documents)
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1, ], times = x[2, ])),
                      use.names = FALSE)
  wcounts <- list(Group.1 = sort(unique(wcountvec)))
  V <- length(wcounts$Group.1)
  if (!isTRUE(all.equal(wcounts$Group.1, 1:V))) {
    stop("Word indices must be sequential integers starting with 1 (every vocabulary term must appear in the corpus).")
  }
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)
  if (length(vocab) != V) stop("Vocab length does not match observed word indices")
  if (!(length(K) == 1 && K == round(K) && K > 1)) {
    stop("K must be a positive integer greater than 1.")
  }

  ## ---- prevalence design matrix (as stm's makeTopMatrix for ~X) ----
  xmat <- cbind(1, as.matrix(X))
  if (anyNA(xmat)) stop("Missing values in prevalence covariates.")
  propSparse <- 1 - sum(xmat != 0) / length(xmat)
  if (!(propSparse < .5 | ncol(xmat) < 50)) {
    xmat <- Matrix::Matrix(xmat, sparse = TRUE)
  }

  ## ---- content covariate ----
  if (!is.null(content)) {
    yvar <- as.factor(content)
    if (any(is.na(yvar))) stop("The content covariate contains missing values.")
    betaindex <- as.numeric(yvar)
  } else {
    betaindex <- rep(1, N)
  }
  A <- length(unique(betaindex))
  if (length(betaindex) != N || nrow(xmat) != N) {
    stop("number of observations in the covariates and documents are not all equal.")
  }

  ## stm::stm() (called without a seed) draws a seed from the current RNG
  ## stream and resets the generator with it.  This is reproduced so that
  ## results obtained with set.seed() are identical to those of sts <= 1.5.
  seed <- floor(stats::runif(1) * 1e7)
  set.seed(seed)

  ## ---- spectral initialization ----
  maxV <- if (V > 10000) 10000 else NULL
  beta <- spectralInit(documents, K, V, maxV = maxV)
  beta <- list(beta = rep(list(beta), A))
  mu <- list(mu = matrix(0, nrow = (K - 1), ncol = 1))
  sigma <- diag(20, nrow = (K - 1))
  lambda <- matrix(0, nrow = N, ncol = (K - 1))

  ## ---- EM (skipped when max.em.its = 0) ----
  if (max.em.its > 0) {
    interactions <- !is.null(content)
    convergence <- NULL
    stopits <- FALSE
    while (!stopits) {
      suffstats <- stmEstep(documents = documents, beta.index = betaindex,
                            update.mu = (!is.null(mu$gamma)),
                            beta = beta$beta, lambda.old = lambda,
                            mu = mu$mu, sigma = sigma)
      sigma.ss <- suffstats$sigma
      lambda <- suffstats$lambda
      beta.ss <- suffstats$beta
      bound.ss <- suffstats$bound
      mu <- opt.mu(lambda = lambda, covar = xmat, enet = 1, ic.k = 2, maxits = 1000)
      sigma <- opt.sigma(nu = sigma.ss, lambda = lambda, mu = mu$mu, sigprior = 0)
      beta <- stmMnreg(beta.ss, K = K, A = A, wcounts = wcounts$x,
                       interactions = interactions)
      convergence <- stmConvergence(bound.ss, convergence, emtol = emtol,
                                    maxits = max.em.its)
      stopits <- convergence$stopits
    }
  }

  logbeta <- lapply(beta$beta, safelog)
  list(eta = lambda, invsigma = solve(sigma),
       beta = list(logbeta = logbeta, kappa = beta$kappa),
       vocab = vocab, settings = list(dim = list(K = K, A = A, V = V, N = N,
                                                 wcounts = wcounts)))
}

## Spectral (anchor word) initialization of the K x V topic-word matrix
## (stm.init with mode = "Spectral").
spectralInit <- function(documents, K, V, maxV = NULL) {
  if (K >= V) stop("Spectral initialization cannot be used for the overcomplete case (K greater than or equal to number of words in vocab)")
  docs <- doc.to.ijv(documents)
  mat <- Matrix::sparseMatrix(docs$i, docs$j, x = docs$v)
  rm(docs)
  wprob <- Matrix::colSums(mat)
  wprob <- wprob / sum(wprob)
  keep <- NULL
  if (!is.null(maxV)) {
    keep <- order(wprob, decreasing = TRUE)[1:maxV]
    mat <- mat[, keep]
    wprob <- wprob[keep]
  }
  Q <- gram(mat)
  Qsums <- rowSums(Q)
  if (any(Qsums == 0)) {
    temp.remove <- which(Qsums == 0)
    if (is.null(keep)) {
      keep <- which(Qsums != 0)
    } else {
      keep <- keep[which(Qsums != 0)]
    }
    Q <- Q[-temp.remove, -temp.remove]
    Qsums <- Qsums[-temp.remove]
    wprob <- wprob[-temp.remove]
  }
  Q <- Q / Qsums
  anchor <- fastAnchor(Q, K = K)
  beta <- recoverL2(Q, anchor, wprob)$A
  if (!is.null(keep)) {
    beta.new <- matrix(0, nrow = K, ncol = V)
    beta.new[, keep] <- beta
    beta.new <- beta.new + .001 / V
    beta <- beta.new / rowSums(beta.new)
  }
  beta
}

## Gram matrix of word co-occurrences
gram <- function(mat) {
  nd <- Matrix::rowSums(mat)
  mat <- mat[nd >= 2, ]
  nd <- nd[nd >= 2]
  divisor <- nd * (nd - 1)
  Q <- Matrix::crossprod(mat / sqrt(divisor)) -
    Matrix::Diagonal(x = Matrix::colSums(mat / divisor))
  as.matrix(Q)
}

## Greedy anchor word selection (stabilized Gram-Schmidt)
fastAnchor <- function(Qbar, K) {
  basis <- c()
  rowSquaredSums <- rowSums(Qbar^2)
  for (i in 1:K) {
    basis[i] <- which.max(rowSquaredSums)
    maxval <- rowSquaredSums[basis[i]]
    normalizer <- 1 / sqrt(maxval)
    Qbar[basis[i], ] <- Qbar[basis[i], ] * normalizer
    innerproducts <- Qbar %*% Qbar[basis[i], ]
    project <- as.numeric(innerproducts) %o% Qbar[basis[i], ]
    project[basis, ] <- 0
    Qbar <- Qbar - project
    rowSquaredSums <- rowSums(Qbar^2)
    rowSquaredSums[basis] <- 0
  }
  basis
}

## RecoverL2 procedure of Arora et al. (exponentiated gradient solver)
recoverL2 <- function(Qbar, anchor, wprob) {
  X <- Qbar[anchor, ]
  XtX <- tcrossprod(X)
  condprob <- vector(mode = "list", length = nrow(Qbar))
  for (i in 1:nrow(Qbar)) {
    if (i %in% anchor) {
      vec <- rep(0, nrow(XtX))
      vec[match(i, anchor)] <- 1
      condprob[[i]] <- vec
    } else {
      y <- Qbar[i, ]
      solution <- expgrad(X, y, XtX)$par
      if (any(solution <= 0)) {
        solution[solution <= 0] <- .Machine$double.eps
      }
      condprob[[i]] <- solution
    }
  }
  weights <- do.call(rbind, condprob)
  A <- weights * wprob
  A <- t(A) / colSums(A)
  list(A = A)
}

## Exponentiated gradient with L2 loss
expgrad <- function(X, y, XtX = NULL, alpha = NULL, tol = 1e-7, max.its = 500) {
  if (is.null(alpha)) alpha <- 1 / nrow(X)
  alpha <- matrix(alpha, nrow = 1, ncol = nrow(X))
  if (is.null(XtX)) XtX <- tcrossprod(X)
  ytX <- y %*% t(X)
  converged <- FALSE
  eta <- 50
  sse.old <- Inf
  its <- 1
  while (!converged) {
    grad <- (ytX - alpha %*% XtX)
    sse <- sum(grad^2)
    grad <- 2 * eta * grad
    maxderiv <- max(grad)
    alpha <- alpha * exp(grad - maxderiv)
    alpha <- alpha / sum(alpha)
    converged <- abs(sqrt(sse.old) - sqrt(sse)) < tol
    if (its == max.its) break
    sse.old <- sse
    its <- its + 1
  }
  list(par = as.numeric(alpha), its = its, converged = converged)
}

## E-step of the STM (Laplace approximation for each document)
stmEstep <- function(documents, beta.index, update.mu, beta, lambda.old, mu,
                     sigma) {
  V <- ncol(beta[[1]])
  K <- nrow(beta[[1]])
  N <- length(documents)
  A <- length(beta)
  if (!update.mu) mu.i <- as.numeric(mu)
  sigma.ss <- diag(0, nrow = (K - 1))
  beta.ss <- vector(mode = "list", length = A)
  for (i in 1:A) {
    beta.ss[[i]] <- matrix(0, nrow = K, ncol = V)
  }
  bound <- vector(length = N)
  lambda <- vector("list", length = N)
  sigobj <- try(chol.default(sigma), silent = TRUE)
  if (inherits(sigobj, "try-error")) {
    sigmaentropy <- (.5 * determinant(sigma, logarithm = TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  for (i in 1:N) {
    doc <- documents[[i]]
    words <- doc[1, ]
    aspect <- beta.index[i]
    init <- lambda.old[i, ]
    if (update.mu) mu.i <- mu[, i]
    beta.i <- beta[[aspect]][, words, drop = FALSE]
    doc.ct <- doc[2, ]
    optim.out <- stats::optim(par = init, fn = stm_lhoodcpp, gr = stm_gradcpp,
                              method = "BFGS", control = list(maxit = 500),
                              doc_ct = doc.ct, mu = mu.i, siginv = siginv,
                              beta = beta.i)
    doc.results <- stm_hpbcpp(optim.out$par, doc_ct = doc.ct, mu = mu.i,
                              siginv = siginv, beta = beta.i,
                              sigmaentropy = sigmaentropy)
    sigma.ss <- sigma.ss + doc.results$eta$nu
    beta.ss[[aspect]][, words] <- doc.results$phis + beta.ss[[aspect]][, words]
    bound[i] <- doc.results$bound
    lambda[[i]] <- c(doc.results$eta$lambda)
  }
  lambda <- do.call(rbind, lambda)
  list(sigma = sigma.ss, beta = beta.ss, bound = bound, lambda = lambda)
}

## M-step for the topic-word distributions with content covariates:
## distributed Poisson (multinomial) regression with an L1 penalty,
## selecting the penalty by AIC (stm's mnreg with default settings).
stmMnreg <- function(beta.ss, K, A, wcounts, interactions = TRUE,
                     alpha = 1, maxit = 1e4, nlambda = 250,
                     lambda.min.ratio = .001, ic.k = 2, thresh = 1e-5) {
  counts <- do.call(rbind, beta.ss)
  if (A == 1) {
    covar <- diag(1, nrow = K)
  } else {
    veci <- 1:nrow(counts)
    vecj <- rep(1:K, A)
    veci <- c(veci, 1:nrow(counts))
    vecj <- c(vecj, rep((K + 1):(K + A), each = K))
    if (interactions) {
      veci <- c(veci, 1:nrow(counts))
      vecj <- c(vecj, (K + A + 1):(K + A + nrow(counts)))
    }
    vecv <- rep(1, length(veci))
    covar <- Matrix::sparseMatrix(veci, vecj, x = vecv)
  }
  m <- log(wcounts) - log(sum(wcounts))
  mult.nobs <- rowSums(counts)
  offset <- log(mult.nobs)
  counts <- split(counts, col(counts))

  subM <- function(x, p) {
    ind <- (x@p[p] + 1):x@p[p + 1]
    rn <- x@i[ind] + 1
    y <- x@x[ind]
    out <- rep(0, length = nrow(x))
    out[rn] <- y
    out
  }

  out <- vector(mode = "list", length = length(counts))
  for (i in 1:length(counts)) {
    offset2 <- m[i] + offset
    mod <- NULL
    while (is.null(mod)) {
      mod <- tryCatch(glmnet::glmnet(x = covar, y = counts[[i]], family = "poisson",
                                     offset = offset2, standardize = FALSE,
                                     intercept = FALSE,
                                     lambda.min.ratio = lambda.min.ratio,
                                     nlambda = nlambda, alpha = alpha,
                                     maxit = maxit, thresh = thresh),
                      warning = function(w) return(NULL),
                      error = function(e) stop(e))
      if (is.null(mod)) nlambda <- nlambda + floor(.2 * nlambda)
    }
    dev <- (1 - mod$dev.ratio) * mod$nulldev
    ic <- dev + ic.k * mod$df
    lambda <- which.min(ic)
    out[[i]] <- subM(mod$beta, lambda)
  }
  coef <- do.call(cbind, out)
  kappa <- split(coef, row(coef))
  linpred <- as.matrix(covar %*% coef)
  linpred <- sweep(linpred, 2, STATS = m, FUN = "+")
  explinpred <- exp(linpred)
  beta <- explinpred / rowSums(explinpred)
  beta <- split(beta, rep(1:A, each = K))
  beta <- lapply(beta, matrix, nrow = K)
  list(beta = beta, kappa = list(m = m, params = kappa))
}

## Convergence check of the STM EM algorithm
stmConvergence <- function(bound.ss, convergence, emtol, maxits) {
  if (is.null(convergence)) {
    convergence <- list(bound = c(), its = 1, converged = FALSE, stopits = FALSE)
  }
  convergence$bound[convergence$its] <- sum(bound.ss)
  if (convergence$its > 1) {
    old <- convergence$bound[convergence$its - 1]
    new <- convergence$bound[convergence$its]
    convergence.check <- (new - old) / abs(old)
    if (emtol != 0 && convergence.check < emtol) {
      convergence$converged <- TRUE
      convergence$stopits <- TRUE
      return(convergence)
    }
  }
  if (convergence$its == maxits) {
    convergence$stopits <- TRUE
    return(convergence)
  }
  convergence$its <- convergence$its + 1
  convergence
}

## Convert documents to triplet (i, j, v) format
doc.to.ijv <- function(documents, fixzeroindex = TRUE) {
  indices <- unlist(lapply(documents, '[', 1, ))
  if ((0 %in% indices) & fixzeroindex) indices <- indices + 1
  counts <- lapply(documents, '[', 2, )
  VsubD <- unlist(lapply(counts, length))
  rowsums <- unlist(lapply(counts, sum))
  docids <- rep(1:length(documents), times = VsubD)
  counts <- unlist(counts)
  list(i = as.integer(docids), j = as.integer(indices), v = as.integer(counts),
       rowsums = as.integer(rowsums))
}

## Convert triplet format to the list-of-documents format
ijv.to.doc <- function(i, j, v) {
  index <- split(j, i)
  index <- lapply(index, as.integer)
  count <- split(v, i)
  count <- lapply(count, as.integer)
  mapply(rbind, index, count, SIMPLIFY = FALSE)
}
