#' Regression Table Estimation
#' 
#' Estimates regression tables for prevalence and sentiment/discourse.
#' 
#' Estimate Gamma coefficients (along with standard errors, p-values, etc.) to  
#' assess how document-level meta-data determine prevalence and sentiment/discourse  
#'  
#' @param alpha the estimated alpha variable values for each document
#' @param mu the mean (fitted) values for alpha based on document-level variables * estimated 
#' Gamma for each document
#' @param sigma the estimated covariance matrix for the alpha parameters
#' @param X the covariates used to estimate the STS model
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
#' X <- model.matrix(~1+out$meta$treatment + out$meta$pid_rep + 
#' out$meta$treatment * out$meta$pid_rep)[,-1]
#' X_seed <- as.matrix(out$meta$treatment)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(X, X_seed, out, numTopics = 3, verbose = FALSE, 
#' parallelize = FALSE, maxIter = 3, initialization = 'anchor')
#' regn_tables <- estimateRegnTables(sts_estimate$alpha, mu = sts_estimate$mu, 
#' sigma = sts_estimate$sigma, X = X)
#' printRegnTables(regn_tables)
#' }
#' @export
estimateRegnTables = function(alpha, mu, sigma, X) {
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
    colnames(coefficients) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    tables[[k]] <- coefficients
  }
  tables
}
