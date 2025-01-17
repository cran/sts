#' Variational EM for the Structural Topic and Sentiment-Discourse (STS) Model
#' 
#' 
#' Estimation of the STS Model using variational EM.
#' The function takes sparse representation of a document-term matrix, covariates 
#' for each document, and an integer number of topics and returns fitted model 
#' parameters. See an overview of functions in the package here: 
#' \code{\link{sts-package}}
#' 
#' This is the main function for estimating the Structural Topic and 
#' Sentiment-Discourse (STS) Model. Users provide a corpus of documents and a 
#' number of topics.  Each word in a document comes from exactly one topic and 
#' each document is represented by the proportion of its words that come from 
#' each of the topics. The document-specific content covariates affect how much 
#' (prevalence) and the way in which a topic is discussed (sentiment-discourse). 
#' 
#' @param prevalence_sentiment A formula object with no response variable or a 
#' design matrix with the covariates. The variables must be 
#' contained in corpus$meta.
#' @param initializationVar A formula with a single variable for use in the initialization of latent sentiment. This argument  
#' is usually the key experimental variable (e.g., review rating binary indicator of experiment/control group).
#' @param corpus The document term matrix to be modeled in a sparse term count matrix with one row
#' per document and one column per term. The object must be a list of with each element 
#' corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is the same format in the \code{\link[stm]{stm}} package. 
#' @param K A positive integer (of size 2 or greater) representing
#' the desired number of topics. 
#' @param maxIter A positive integer representing the max number of VEM iterations allowed.
#' @param convTol Convergence tolerance for the variational EM estimation algorithm; Default value = 1e-5.
#' @param initialization Character argument that allows the user to specify an initialization
#' method. The default choice, \code{"anchor"} to initialize prevalence according to anchor words and 
#' the key experimental covariate identified in argument \code{initializationVar}. One can also use 
#' \code{"stm"}, which uses a fitted STM model (Roberts et al. 2014, 2016) 
#' to initialize coefficients related to prevalence and sentiment-discourse. 
#' @param kappaEstimation A character input specifying how kappa should be estimated. \code{"lasso"} allows for 
#' penalties on the L1 norm.  We estimate a regularization path and then select the optimal
#' shrinkage parameter using AIC. \code{"adjusted"} (default) utilizes the lasso penalty with an adjusted aggregated Poisson regression. 
#' All options use an approximation framework developed in Taddy (2013) called
#' Distributed Multinomial Regression which utilizes a factorized poisson
#' approximation to the multinomial.  See Li and Mankad (2024) on the implementation here.  
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.  
#' @param parallelize A logical flag indicating whether to parallelize the estimation using all but one CPU cores on your local machine. 
#' @param stmSeed A prefit STM model object to initialize the STS model. Note this is ignored unless initialization = "stm"
#' @return An object of class sts 
#' 
#' \item{alpha}{Estimated prevalence and sentiment-discourse values for each document and topic} 
#' \item{gamma}{Estimated regression coefficients that determine prevalence and sentiment/discourse for each topic} 
#' \item{kappa}{Estimated kappa coefficients that determine sentiment-discourse and the topic-word distributions} 
#' \item{sigma_inv}{Inverse of the covariance matrix for the alpha parameters}
#' \item{sigma}{Covariance matrix for the alpha parameters} 
#' \item{elbo}{the ELBO at each iteration of the estimation algorithm}
#' \item{mv}{the baseline log-transformed occurrence rate of each word in the corpus}
#' \item{runtime}{Time elapsed in seconds} 
#' \item{vocab}{Vocabulary vector used} 
#' \item{mu}{Mean (fitted) values for alpha based on document-level variables * estimated 
#' Gamma for each document}
#' 
#' 
#' @seealso  \code{\link{estimateRegns}}
#' @references 
#' Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013)
#' "The structural topic model and applied social science." In Advances in
#' Neural Information Processing Systems Workshop on Topic Models: Computation,
#' Application, and Evaluation. 
#' 
#' Roberts M., Stewart, B. and Airoldi, E. (2016) "A model of text for
#' experimentation in the social sciences" Journal of the American Statistical
#' Association.
#' 
#' Chen L. and Mankad, S. (2024) "A Structural Topic and Sentiment-Discourse Model
#' for Text Analysis" Management Science.
#' @examples
#' #An example using the Gadarian data from the stm package.  From Raw text to 
#' # fitted model using textProcessor() which leverages the tm Package
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 1, verbose = FALSE)
#' @export
sts = function(prevalence_sentiment, initializationVar, corpus, K, maxIter = 100, convTol = 1e-5, initialization = "anchor", kappaEstimation = "adjusted", verbose = TRUE, parallelize = FALSE, stmSeed = NULL) {

  # if (!is.null(seed)) {set.seed(seed)}
  if(inherits(prevalence_sentiment,"formula")) {
    X <- model.matrix(prevalence_sentiment, data = corpus$meta)[,-1, drop = FALSE]
  } else {
    X <- prevalence_sentiment
  }
  if(inherits(initializationVar,"formula")) {
    X_seed <- model.matrix(initializationVar, data = corpus$meta)[,-1, drop = FALSE]
  } else {
    X_seed <- initializationVar
  }
  
  if (parallelize) {
    # Register the parallel backend
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl, varlist = c("lpbdcpp", "lgaecpp", "esthcpp"))
    
    # Function to handle each iteration of the loop
    process_document <- function(d) {
      cdv_ind <- corpus$documents[[d]][1,]
      cdv <- corpus$documents[[d]][2,]
      document <- corpus$documents[[d]]
      
      infer <- optim(par = alpha.est[d,], fn = lpbdcpp, gr = lgaecpp,
                     method = "BFGS", control = list(fnscale = -1, maxit=500), hessian = FALSE,
                     Sigma_Inv=Sigma_Inv.est,
                     kappa_t=kappa.est$kappa_t, kappa_s=kappa.est$kappa_s, mu_d = mu$mu[,d], doc = document, V=V, mv = mv)
      
      hessian <- esthcpp(alpha_d = infer$par, kappa_t=kappa.est$kappa_t, kappa_s=kappa.est$kappa_s, Sigma_Inv = Sigma_Inv.est, doc = document, V=V, mv = mv)
      h.chol <- tryCatch(expr = chol(hessian), error = function(x) return(NULL))
      if (is.null(h.chol)) {
        dvec <- diag(hessian)
        magnitudes <- rowSums(abs(hessian)) - abs(dvec)
        dvec <- mapply(max, dvec, magnitudes)
        diag(hessian) <- dvec
        h.chol <- chol(hessian)
      }
      sigma.ss <- chol2inv(h.chol)
      elbo <- infer$value + 0.5*(sum(log(diag(h.chol))) - log(det(sigma)))
      
      # if(verbose && d%%ctevery==0) cat(".")
      
      list(sigma.ss = sigma.ss, elbo = elbo, infer = infer)
    }
    
    
    
  }

  conv_criteria <- convTol
  topicreportevery <- 10
  estimation = kappaEstimation
  nc_X <- ncol(X)
  D <- length(corpus$documents)
  V <- length(corpus$vocab)
  c_d_bar <- unlist(lapply(corpus$documents, function(x) sum(x[2,])))
  wcountvec <- unlist(lapply(corpus$documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
  wcounts <- table(wcountvec)
  mv <- as.numeric(log(wcounts/sum(wcounts)))
  
  K <- K
  if (verbose) {cat("\nStarting Initialization...")}
  
  
  if (verbose) {cat("Beta...")}
  if (initialization == "anchor") {
    ## anchor words initialization
    mod.out <- suppressWarnings(stm(corpus$documents, corpus$vocab, K, verbose = FALSE,
                                    prevalence=~X, data=data.frame(X), max.em.its = 0, LDAbeta = FALSE))
  } else if (initialization == "stm") {
    ## stm initialization
    mod.out <- stmSeed
    if (is.null(stmSeed))
      mod.out <- suppressWarnings(stm(corpus$documents, corpus$vocab, K, verbose = FALSE, content = X_seed,
                                      prevalence=~X, data=data.frame(X), max.em.its = 100, LDAbeta = FALSE))
  }
  
  if (verbose) {cat("Sigma, mu, and alpha...")}
  ## define the other parameters
  Sigma_Inv.est <- diag(1/20, 2*K-1)
  Sigma_Inv.est[1:(K-1),1:(K-1)] <- mod.out$invsigma
  sigma <- solve(Sigma_Inv.est)
  alpha.est <- matrix(0, D, 2*K-1)
  alpha.est[,1:(K-1)] <- mod.out$eta
  alpha.est[,1:K+K-1] <- matrix(X_seed - mean(X_seed), byrow = FALSE)
  # alpha.est[,1:K+K-1] <- matrix(1, byrow = FALSE)
  mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000)
  Gamma.est <- mu$gamma
  
  ## start sequence to get kappa
  if (verbose) {cat("Kappa...")}
  rho <- 0    
  kappa.est <- list(kappa_t = matrix(0, V, K), kappa_s = matrix(0, V, K))
  
  ## create groups for aggregation (to eventually estimate kappa)
  X_seed_orig <- X_seed
  group <- rep(1, D)
  vals <- sort(unique(X_seed))
  numGroups <- length(vals)
  for (g in 1:numGroups) {
    group[X_seed == vals[g]] <- g
  } 
  
  ## calculate phi
  phi <- vector("list", length = numGroups)
  phiD <- matrix(0, nrow=D, ncol = K)
  for (g in 1:numGroups) {
    phi[[g]] <- matrix(data = 0, nrow = V, ncol = K)
  }  
  BZ <- t(exp(mod.out$beta$logbeta[[1]]))
  BZ <- apply(BZ, 2, function(x) x / rowSums(BZ))
  for (d in 1:D) {
    cdv_ind <- corpus$documents[[d]][1,]
    cdv <- corpus$documents[[d]][2,]
    phi[[group[d]]][cdv_ind,] <- phi[[group[d]]][cdv_ind,] + cdv*BZ[cdv_ind,]
    phiD[d,] <- phiD[d,] + colSums(cdv*BZ[cdv_ind,, drop = FALSE])
  }

  ## estimate kappa 
  if (initialization == "anchor") {
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, c_d_bar = c_d_bar, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, phiD =phiD, parallelize = parallelize)
    kappa.est$kappa_t[which(kappa.est$kappa_t == 0)] <- sign(rnorm(sum(kappa.est$kappa_t == 0)))*1e-4
    kappa.est$kappa_s[which(kappa.est$kappa_s == 0)] <- sign(rnorm(sum(kappa.est$kappa_s == 0)))*1e-4
  } else {
    kappa.est$kappa_t <- do.call(cbind,mod.out$beta$kappa$params[1:K])
    
    kappa_content <- do.call(cbind, mod.out$beta$kappa$params[(K+1):(K+numGroups+1)])
    
    tindx <- rep(1:K, each = numGroups)
    kappa_interaction <- do.call(cbind, mod.out$beta$kappa$params[tindx])
    
    kappa_s_neg <- kappa_interaction[,seq(1,K*numGroups, by = numGroups)]
    kappa_s_pos <- kappa_interaction[,seq(numGroups,K*numGroups, by = numGroups)]
    #mod.out$beta$kappa$params[[(K+1)]] +
    #do.call(cbind, mod.out$beta$kappa$params[(K+numGroups+1):(K+numGroups+K)])
    
    # print((K+numGroups+1):(K+numGroups+K))  
    # A <- model$settings$dim$A ## numGroups
    # anames <- model$settings$covariates$yvarlevels
    # i1 <- K + 1
    # i2 <- K + A
    # intnums <- (i2 + 1):nrow(labs)
    # out$topics <- labs[topics, , drop = FALSE]
    # out$covariate <- labs[i1:i2, , drop = FALSE]
    # rownames(out$covariate) <- anames
    # if (model$settings$kappa$interactions) {
    #   tindx <- rep(1:K, each = A)
    #   intnums <- intnums[tindx %in% topics]
    #   out$interaction <- labs[intnums, , drop = FALSE]
    # }
    
    # print(K)
    # print(numGroups)
    
    # kappa_s_pos <- mod.out$beta$kappa$params[[(K+numGroups)]] +
    # do.call(cbind, mod.out$beta$kappa$params[(K+numGroups+K+1):length(mod.out$beta$kappa$params)])
    
    # print((K+numGroups+K+1):length(mod.out$beta$kappa$params))  
    
    # print(dim(kappa_s_pos))
    # print(dim(kappa_s_neg))
    
    kappa.est$kappa_s <- kappa_s_pos - kappa_s_neg
    
  }
  
  elbo_vec <- achange_vec <- kchange_vec <- rep(NA, maxIter)
  ntokens <- sum(mod.out$settings$dim$wcounts$x)
  converge <- FALSE
  iter <- 1
  if (verbose) {cat("Done.\nStarting Estimation...\n")}
  timer_full <- proc.time()[3]  

  while (iter <= maxIter & !converge) {
    t1 <- proc.time()[3]
    alpha.est.old <- alpha.est ## DEBUG
    kappa.est.old <- kappa.est ## DEBUG
    sigma.ss <- diag(0, nrow=2*K-1)
    phi <- vector("list", length = numGroups)
    phiD <- matrix(0, nrow=D, ncol = K)
    for (g in 1:numGroups) {
      phi[[g]] <- matrix(data = 0, nrow = V, ncol = K)
    }  
    elbo <- 0      
    BZ <- matrix(0, V, K)
    ctevery <- ifelse(D>100, floor(D/100), 1)
    
    if (!parallelize) {
      for (d in 1:D) {
        ###################################################### THIS CODE CAN BE OPTIMIZED IN C++ ######################################################
        cdv_ind <- corpus$documents[[d]][1,]
        cdv <- corpus$documents[[d]][2,]
        document <- corpus$documents[[d]]
        # if (cpp) {
          infer <- optim(par = alpha.est[d,], fn = lpbdcpp, gr = lgaecpp,
                         method = "BFGS", control = list(fnscale = -1, maxit=500), hessian = FALSE,
                         Sigma_Inv=Sigma_Inv.est,
                         kappa_t=kappa.est$kappa_t, kappa_s=kappa.est$kappa_s, mu_d = mu$mu[,d], doc = document, V=V, mv = mv
          )
        # }
        # else {
        #   infer <- optim(par = alpha.est[d,], fn = log_posterior_byDoc.optim, gr = lapl_grad_alpha_eta.optim,
        #                  method = "BFGS", control = list(fnscale = -1, maxit=500), hessian = FALSE,
        #                  Gamma=Gamma.est, Sigma_Inv=Sigma_Inv.est,
        #                  kappa=kappa.est, mu_d = mu$mu[,d], doc = document, V=V, mv = mv
        #   )
        # }
        alpha.est[d,] <- infer$par
        # if (cpp) {
          hessian <- esthcpp(alpha_d = infer$par, kappa_t=kappa.est$kappa_t, kappa_s=kappa.est$kappa_s, Sigma_Inv = Sigma_Inv.est, doc = document, V=V, mv = mv)
        # }
        # else{
          # hessian <- estimateHessian(alpha.d = infer$par, kappa = kappa.est, Sigma_Inv = Sigma_Inv.est, doc = document, V=V, mv = mv)
        # }
        h.chol <-  tryCatch(expr = chol(hessian), error = function(x) return(NULL))
        if(is.null(h.chol)){
          dvec <- diag(hessian)
          magnitudes <- rowSums(abs(hessian)) - abs(dvec)
          dvec <- mapply(max, dvec, magnitudes)
          diag(hessian) <- dvec
          h.chol <-  chol(hessian)
        }
        sigma.ss <- sigma.ss + chol2inv(h.chol)
        elbo <- elbo + infer$value + 0.5*(sum(log(diag(h.chol))) - log(det(sigma)))
        ###################################################### STOP ######################################################
  
        expeta <- exp(c(infer$par[1:(K-1)],0))
        theta <- expeta/sum(expeta)
        for (k in 1:K) {
          BZ[,k] <- theta[k]*softmax(mv + kappa.est$kappa_t[,k] + kappa.est$kappa_s[,k] * alpha.est[d,K+k-1])
        }
        BZ <- apply(BZ, 2, function(x) x / rowSums(BZ))
        phi[[group[d]]][cdv_ind,] <- phi[[group[d]]][cdv_ind,] + cdv*BZ[cdv_ind,]
        phiD[d,] <- phiD[d,] + colSums(cdv*BZ[cdv_ind,, drop = FALSE])

        if(verbose && d%%ctevery==0) cat(".")
      }
  
      # print(summary(diff))
      if(verbose) {
        if(verbose) cat("\n") #add a line break for the next message.
        msg <- sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
        cat(msg)
      }
    } else {
      if(verbose) {
        cat("Starting E-step...")
      }
      
      
      
      
      # Use foreach for parallel processing
      i = 1
      results <- foreach(i = 1:D) %dopar% {
        # dyn.load("STS")
        process_document(i)
      }
      
      
      # Parallelize using mclapply
      # results <- mclapply(1:D, process_document, mc.cores = detectCores() - 1)
      
      # Create a cluster using available cores
      # cl <- makeCluster(detectCores() - 1)
      # clusterExport(cl, varlist = c("out", "alpha.est", "lpbdcpp", "lgaecpp", "Sigma_Inv.est", "kappa.est", "mu", "V", "mv", "process_document", "esthcpp", "sigma"))
      # results <- parLapply(cl, 1:D, process_document)
      # stopCluster(cl)
      
      # Aggregate results
      for (d in 1:D) {
        sigma.ss <- sigma.ss + results[[d]]$sigma.ss
        elbo <- elbo + results[[d]]$elbo
        infer <- results[[d]]$infer
        alpha.est[d,] <- infer$par
        
        cdv_ind <- corpus$documents[[d]][1,]
        cdv <- corpus$documents[[d]][2,]
        expeta <- exp(c(infer$par[1:(K-1)], 0))
        theta <- expeta / sum(expeta)
        for (k in 1:K) {
          BZ[,k] <- theta[k] * softmax(mv + kappa.est$kappa_t[,k] + kappa.est$kappa_s[,k] * infer$par[K+k-1])
        }
        BZ <- apply(BZ, 2, function(x) x / rowSums(BZ))
        phi[[group[d]]][cdv_ind,] <- phi[[group[d]]][cdv_ind,] + cdv * BZ[cdv_ind,]
        phiD[d,] <- phiD[d,] + colSums(cdv*BZ[cdv_ind,, drop = FALSE])
      }
      if(verbose) {
        # if(verbose) cat("\n") #add a line break for the next message.
        msg <- sprintf("Completed E-Step (%d seconds). \n", floor((proc.time()-t1)[3]))
        cat(msg)
      }
    }
    t1 <- proc.time()
    
    
    # Now do the maximum likelihood updates 
    mu <- opt.mu(lambda=alpha.est, covar=cbind(1,X), enet=0, ic.k=NULL, maxits=1000)
    Gamma.est <- mu$gamma
    
    sigma <- opt.sigma(nu=sigma.ss, lambda=alpha.est,
                       mu=mu$mu, sigprior=0)
    sigma.chol <- chol(sigma)
    Sigma_Inv.est <- chol2inv(sigma.chol)
    
    kappa.est <- opt.kappa(phi = phi, alphaS = alpha.est[,1:K+K-1], kappa = kappa.est, c_d_bar = c_d_bar, estimation = estimation, numGroups = numGroups, group = group, V = V, mv = mv, phiD =phiD, parallelize = parallelize)
    
    if(verbose) {
      #M-step message
      timer <- floor((proc.time()-t1)[3])
      msg <- ifelse(timer>1,
                    sprintf("Completed M-Step (%d seconds). \n", floor((proc.time()-t1)[3])),
                    "Completed M-Step. \n")
      
      cat(msg)
    }
    
    elbo_vec[iter] <- elbo
    if (iter == 1) {
      msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f) \n", 
                     iter, elbo/ntokens) 
      if(verbose) {
        cat(msg)
      }
    } else {
      new <- elbo_vec[iter]
      old <- elbo_vec[iter-1]
      delta <- (new - old)/(0.1 + abs(old))
      converge <- ifelse(delta < conv_criteria, TRUE, FALSE)
      
      achange <- mean(abs(alpha.est - alpha.est.old)/(0.1 + abs(alpha.est.old)))
      kchange <- mean(abs(c(kappa.est$kappa_t, kappa.est$kappa_s) - c(kappa.est.old$kappa_t, kappa.est.old$kappa_s))/(0.1 + abs(c(kappa.est.old$kappa_t, kappa.est.old$kappa_s))))
      
      achange_vec[iter] <- achange
      kchange_vec[iter] <- kchange
      
      if(verbose) {
        msg <- sprintf("Completing Iteration %i (approx. per word bound = %.3f, relative change = %.3e;\n\t\t\tavg rel change in alpha = %.3e; avg rel change in kappa = %.3e) \n",
                       iter, new/ntokens, delta, achange, kchange)
        cat(msg)
      }
      if (verbose && iter %% topicreportevery==0) {
        printTopWords(list(alpha = alpha.est, kappa = kappa.est, sigma = sigma, mv = mv, vocab = mod.out$vocab))
      }        
    }     
    iter <- iter + 1
  }

  final.est <- list(alpha = alpha.est, gamma = Gamma.est, kappa = kappa.est, sigma_inv = Sigma_Inv.est, sigma = sigma, elbo = elbo_vec[1:iter], mv = mv, runtime = proc.time()[3] - timer_full, vocab = mod.out$vocab, mu = mu$mu)
  if (parallelize) {
    # Stop the cluster when done
    parallel::stopCluster(cl)
  }
  class(final.est) <- "STS"

  return(final.est)
}