#' Regression Table Estimation
#' 
#' Estimates regression tables for prevalence and sentiment/discourse.
#' 
#' Estimate Gamma coefficients (along with standard errors, p-values, etc.) to  
#' assess how document-level meta-data determine prevalence and sentiment/discourse  
#'  
#' @param object an sts object
#' @param prevalence_sentiment A formula object with no response variable or a 
#' design matrix with the covariates. If a formula, the variables must be 
#' contained in corpus$meta.
#' @param corpus The document term matrix to be modeled in a sparse term count matrix with one row
#' per document and one column per term. The object must be a list of with each element 
#' corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is the same format in the \code{\link[stm]{stm}} package. 
#' 
#' @return a list of tables with regression coefficient estimates. The first 
#' <num-topic> elements pertain to prevalence; the latter  <num-topic> elements 
#' pertain to sentiment-discourse.
#' 
#' @examples \donttest{
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' regns <- estimateRegns(sts_estimate, ~treatment*pid_rep, out)
#' printRegnTables(x = regns)
#' }
#' @export
estimateRegns = function(object, prevalence_sentiment, corpus) {
  
  alpha <- object$alpha
  mu <- object$mu
  sigma <- object$sigma
  
  if(inherits(prevalence_sentiment,"formula")) {
    X <- model.matrix(prevalence_sentiment, data = corpus$meta)[,-1, drop = FALSE]
  } else {
    X <- prevalence_sentiment
  }

  D <- nrow(alpha)
  K <- ceiling(ncol(alpha)/2)

  row.lse <- function(mat) {
    matrixStats::rowLogSumExps(mat)
  }
  
  # A function for performing simple linear regression with a cached QR decomposition
  # this should be lighter weight than lm().  Works with summary.qr.lm() to give
  # vcov calculations etc.
  qr.lm <- function(y, qx) {
    if(length(y)!=nrow(qx$qr)) {
      #probably don't match because of a prior
      if(length(y)!=(nrow(qx$qr)-ncol(qx$qr))) stop("number of covariate observations does not match number of docs")
      #if it passes this check its the prior. thus
      y <- c(y,rep(0, ncol(qx$qr)))
    }
    beta <- solve.qr(qx, y)
    residuals <- qr.resid(qx,y)
    fitted.values <- qr.fitted(qx,y)
    df.residual <- length(fitted.values) - qx$rank
    out <- list(coefficients=beta, residuals=residuals, 
                fitted.values=fitted.values, 
                df.residual=df.residual, rank=qx$rank, qr=qx)
    out 
  }
  #this function rewrites the summary.lm() function
  # to calculate from our reduced regression
  summary.qr.lm <- function (object) {
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    
    Qr <- object$qr
    n <- nrow(Qr$qr)
    p1 <- 1L:p
    r <- z$residuals
    f <- z$fitted.values
    
    mss <- ifelse(attr(z$terms, "intercept"), sum((f - mean(f))^2), sum(f^2)) 
    rss <- sum(r^2)
    
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    sigma <- sqrt(resvar)
    list(est=est, vcov=(sigma^2 * R))
  }
  
  if(ncol(mu)==1) { 
    #if there is only one global mean vector we avoid storing them all, thus calc done with a sweep
    covariance <- crossprod(sweep(alpha, 2, STATS=as.numeric(mu), FUN="-"))
  } else {
    #the typical calculation when there are frequency covariates
    covariance <- crossprod(alpha-t(mu)) 
  }
  #rescale by the number of documents
  covariance <- covariance/nrow(alpha)
  #subtract off the effect from the global covariance
  Sigma <- sigma - covariance
  
  nsims <- 25
  nsims_se <- 250
  tables <- vector(mode="list", length=2*K)
  storage <- vector(mode="list", length=2*K)
  for (i in 1:nsims) {
    out <- vector(mode="list",length=D) 
    for (d in 1:D) {
      mat <- mvtnorm::rmvnorm(1, alpha[d,],sigma,checkSymmetry = FALSE)
      eta <- cbind(mat[,1:(K-1),drop = FALSE],1)
      alpha_s <- mat[,1:K+K-1]
      mat <- c(exp(eta - row.lse(eta)), alpha_s)
      out[[d]] <- mat
    }
    thetasims <- do.call(rbind, out)
    qx <- qr(cbind(1, X))
    for (k in 1:(2*K)){
      # lm.mod <- lm(thetasims[,k]~ qx)
      # storage[[which(k==K)]][[i]] <- list(coef=coef(lm.mod),vcov=vcov(lm.mod))
      lm.mod <- qr.lm(thetasims[,k], qx)
      storage[[k]][[i]] <- summary.qr.lm(lm.mod)
    }
  }
  for (k in 1:(2*K)) {
    sims <- lapply(storage[[k]], function(x) mvtnorm::rmvnorm(nsims_se, x$est, x$vcov))
    sims <- do.call(rbind,sims)
    est<- colMeans(sims)
    se <- sqrt(apply(sims,2, stats::var))
    tval <- est/se
    rdf <- nrow(alpha) - length(est)
    p <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
    
    coefficients <- cbind(est, se, tval, p)
    rownames(coefficients) <- attr(storage[[1]][[1]]$est, "names") 
    rownames(coefficients)[1] <- "Intercept"
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    tables[[k]] <- coefficients
  }
  tables
}
