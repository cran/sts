## Estimates Gamma and fitted values called mu -- code based on the STM packages
# main method up top, regression-implementations below.
opt.mu <- function(lambda, covar=NULL, enet=NULL, ic.k=2,
                   maxits=1000) {

  #Variational Linear Regression with a Gamma hyperprior
    gamma <- vector(mode="list",length=ncol(lambda))
    Xcorr <- crossprod(covar)
    for (i in 1:ncol(lambda)) {
      gamma[[i]] <- vb.variational.reg(Y=lambda[,i], X=covar, Xcorr=Xcorr, maxits=maxits) 
    }
    gamma <- do.call(cbind,gamma)
    mu<- t(covar%*%gamma)
    #if its not a regular matrix,coerce it as it won't be sparse.
    if(!is.matrix(mu)) {
      mu <- as.matrix(mu)
    }
    return(list(mu=mu, gamma=gamma))
}

#Variational Linear Regression with a Half-Cauchy hyperprior 
# (Implementation based off the various LMM examples from Matt Wand)
# This code is intended to be passed a Matrix object
vb.variational.reg <- function(Y,X, b0=1, d0=1, Xcorr=NULL, maxits=1000) {
  if(is.null(Xcorr)) Xcorr <- crossprod(X)
  XYcorr <- crossprod(X,Y) 
  
  an <- (1 + nrow(X))/2
  D <- ncol(X)
  N <- nrow(X)
  w <- rep(0, ncol(X))
  error.prec <- 1 #expectation of the error precision
  converge <- 1000
  cn <- ncol(X) # - 1 for the intercept and +1 in the update cancel
  dn <- 1
  Ea <- cn/dn #expectation of the precision on the weights
  ba <- 1
  
  ct <- 1
  while(converge>.0001) {
    w.old <- w
    
    #add the coefficient prior.  Form depends on whether X is a Matrix object or a regular matrix.
    if(is.matrix(X)) {
      ppmat <- diag(x=c(0, rep(as.numeric(Ea), (D-1))),nrow=D) 
    } else {
      ppmat <- Matrix::Diagonal(n=D, x=c(0, rep(as.numeric(Ea), (D-1))))
    }
    invV <- error.prec*Xcorr + ppmat
    #if its a plain matrix its faster to use the cholesky, otherwise just use solve
    if(is.matrix(invV)) {
      V <- chol2inv(chol(invV))
    } else {
      #Matrix package makes this faster even when its non-sparse
      V <- solve(invV)     
    }
    w <- error.prec*V%*%XYcorr
    
    # parameters of noise model (an remains constant)
    sse <- sum((X %*% w - Y)^ 2)
    bn <- .5*(sse + sum(diag(Xcorr%*%V))) + ba
    error.prec <- an/bn
    ba <- 1/(error.prec + b0)
    
    #subtract off the intercept while working out the hyperparameters
    # for the coefficients
    w0 <- w[1]
    w <- w[-1]
    da <- 2/(Ea + d0)
    dn <- 2*da + (crossprod(w) + sum(diag(V)[-1]))
    Ea <- cn / dn
    #now combine the intercept back in 
    w <- c(w0,w)
    ct <- ct + 1
    if(ct > maxits) {
      stop("Prevalence/Sentiment regression failing to converge within iteration limit.")
    }
    converge <- sum(abs(w-w.old))
  }
  return(w)
}


