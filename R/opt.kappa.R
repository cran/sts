# Kappa Estimation
# 
# Penalized multinomial regression methods for estimating the kappa coefficents using the distributed poisson approach. 
# By default, via the main STS function, estimation is performed with lasso.
opt.kappa <- function(phi, alphaS, kappa, c_d_bar, estimation = estimation, numGroups, group, V, mv, phiD, parallelize = parallelize) {
  low <- -Inf
  upp <- Inf
  K <- ncol(alphaS)
  penalty.factor <- c(rep(1,K), rep(1,K))
  
  phi_sum <- safelog(unlist(lapply(phi, colSums))) # results in numGroups*K vector
  x1_agg <- matrix(rep(diag(K), numGroups), ncol = K, byrow = TRUE)
  
  alpha_s <- lapply(split(alphaS, group), matrix, ncol = K) ## numGroups by topics
  alpha_agg <- Reduce('rbind', lapply(alpha_s, function(m) diag(colMeans(m))))
  
  ## new code ##
  # Create an empty matrix to store the results
  # alpha_agg <- matrix(0, nrow = length(unique(group)) * K, ncol = K)
  # 
  # # Compute colMeans and store the results directly in alpha_agg
  # group_indices <- split(seq_len(nrow(alphaS)), group)
  # row_idx <- 1
  # for (g in seq_along(group_indices)) {
  #   group_matrix <- matrix(alphaS[group_indices[[g]], ], ncol = K)
  #   alpha_agg[row_idx:(row_idx + K - 1), ] <- diag(colMeans(group_matrix))
  #   row_idx <- row_idx + K
  # }
  ## end new code ##
  
  if (!parallelize) {
    for (v in 1:V) {
      y_agg <- unlist(lapply(phi, function(m)m[v,]))
      offset_agg <- phi_sum + mv[v]
      if (estimation == "lasso") {
        mod <- NULL
        inner_count <- 1
        while(is.null(mod)) {
          nlambda <- 250
          mod <- tryCatch(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                         offset=offset_agg, standardize=FALSE,
                                         intercept=FALSE,
                                         lambda.min.ratio=0.001,
                                         nlambda=250, alpha=1, lower.limits = low, upper.limits = upp,
                                         maxit=10000, thresh=1e-05, 
                                         penalty.factor = penalty.factor),
                          warning=function(w) return(NULL),
                          error=function(e) stop(e))
          #if it didn't converge, increase nlambda paths by 50%
          if(is.null(mod)) nlambda <- nlambda + floor(50*nlambda)
          inner_count <-inner_count + 1
          #if it didn't converge after 10, just take the estimates as is
          if (inner_count == 10) {
            mod <- suppressWarnings(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                                   offset=offset_agg, standardize=FALSE,
                                                   intercept=FALSE,
                                                   lambda.min.ratio=0.001,
                                                   nlambda=250, alpha=1, lower.limits = low, upper.limits = upp,
                                                   maxit=10000, thresh=1e-05, 
                                                   penalty.factor = penalty.factor))
          }
        }
        dev <- (1-mod$dev.ratio)*mod$nulldev
        ic <- dev + 2*mod$df ## AIC
        penalty <- which.min(ic)
        coef <- mod$beta[,penalty] #return coefficients
      } else if (estimation == "adjusted") {
        coef <- rep(0, 2*K)
        phiD_list <- lapply(split(phiD, group), matrix, ncol = K) 
        alpha_s <- lapply(split(alphaS, group), matrix, ncol = K) ## numGroups by topics
        G <- length(alpha_s)
        
        # coef <- c(kappa$kappa_t[v,], kappa$kappa_s[v,])
        ## phi is a list of length numGroups, each element is V \times K
        # phiSum_byG <- lapply(phi, function(x) diag(colSums(x))) ## results in list of length numGroups with vectors of length K 
        # delta_converge <- NA
        kk <- 1
        max_inner_iter <- 10
        weighted_alpha <- list()
        while (kk <= max_inner_iter) {
          coef_old <- coef
          
          for (t in 1:G){
            # w[[t]] <- apply(log(phiD_list[[t]]) + (alpha_s[[t]] %*% diag(coef[1:K+K])), 2, softmax)
            weighted_alpha[[t]] <- alpha_s[[t]] * apply(safelog(phiD_list[[t]]) + (alpha_s[[t]] %*% diag(coef[1:K+K])), 2, softmax)
          }
          alpha_agg <- Reduce('rbind', lapply(weighted_alpha, function(x) diag(colSums(x))))

          delta <- coef - coef_old 
          delta_converge <- mean(abs((delta))/abs(0.1+coef_old))
          
          if(delta_converge <= 1e-4) {break}
          kk <- kk + 1
          
        }
        ###### LASSO PORTION ######
        mod <- NULL
        inner_count <- 1
        while(is.null(mod)) {
          nlambda <- 250
          mod <- tryCatch(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                         offset=offset_agg, standardize=FALSE,
                                         intercept=FALSE,
                                         lambda.min.ratio=0.001,
                                         nlambda=nlambda, alpha=1,lower.limits = low, upper.limits = upp,
                                         maxit=10000, thresh=1e-05, 
                                         penalty.factor = penalty.factor),
                          warning=function(w) return(NULL),
                          error=function(e) stop(e))
          #if it didn't converge, increase nlambda paths by 50%
          if(is.null(mod)) nlambda <- nlambda + floor(50*nlambda)
          inner_count <-inner_count + 1
          #if it didn't converge after 10, just take the estimates as is
          if (inner_count == 10) {
            mod <- suppressWarnings(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                                   offset=offset_agg, standardize=FALSE,
                                                   intercept=FALSE,
                                                   lambda.min.ratio=0.001,
                                                   nlambda=nlambda, alpha=1,lower.limits = low, upper.limits = upp,
                                                   maxit=10000, thresh=1e-05, 
                                                   penalty.factor = penalty.factor))
          }
        }
        dev <- (1-mod$dev.ratio)*mod$nulldev
        ic <- dev + 2*mod$df ## AIC
        penalty <- which.min(ic)
        coef <- mod$beta[,penalty] #return coefficients
        ###### END LASSO PORTION ######
      } 
      kappa$kappa_t[v,] <- (coef[1:K] + kappa$kappa_t[v,])/2
      kappa$kappa_s[v,] <- (coef[1:K+K] + kappa$kappa_s[v,])/2
      # kappa$kappa_t[v,] <- coef[1:K]
      # kappa$kappa_s[v,] <- coef[1:K+K]
    }
  } else {
    process_word <- function(v) {
      y_agg <- unlist(lapply(phi, function(m)m[v,]))
      offset_agg <- phi_sum + mv[v]
      if (estimation == "lasso") {
        mod <- NULL
        inner_count <- 1
        while(is.null(mod)) {
          nlambda <- 250
          mod <- tryCatch(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                         offset=offset_agg, standardize=FALSE,
                                         intercept=FALSE,
                                         lambda.min.ratio=0.001,
                                         nlambda=250, alpha=1, lower.limits = low, upper.limits = upp,
                                         maxit=10000, thresh=1e-05, 
                                         penalty.factor = penalty.factor),
                          warning=function(w) return(NULL),
                          error=function(e) stop(e))
          #if it didn't converge, increase nlambda paths by 50%
          if(is.null(mod)) nlambda <- nlambda + floor(50*nlambda)
          inner_count <-inner_count + 1
          #if it didn't converge after 10, just take the estimates as is
          if (inner_count == 10) {
            mod <- suppressWarnings(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                                   offset=offset_agg, standardize=FALSE,
                                                   intercept=FALSE,
                                                   lambda.min.ratio=0.001,
                                                   nlambda=250, alpha=1, lower.limits = low, upper.limits = upp,
                                                   maxit=10000, thresh=1e-05, 
                                                   penalty.factor = penalty.factor))
          }
        }
        dev <- (1-mod$dev.ratio)*mod$nulldev
        ic <- dev + 2*mod$df ## AIC
        penalty <- which.min(ic)
        coef <- mod$beta[,penalty] #return coefficients
      } else if (estimation == "adjusted"){
        coef <- rep(0, 2*K)
        
        # coef <- c(kappa$kappa_t[v,], kappa$kappa_s[v,])
        ## phi is a list of length numGroups, each element is V \times K
        # phiSum_byG <- lapply(phi, function(x) diag(colSums(x))) ## results in list of length numGroups with vectors of length K 
        # delta_converge <- NA
        kk <- 1
        max_inner_iter <- 1
        weighted_alpha <- alpha_s#list()
        while (kk <= max_inner_iter) {
          coef_old <- coef
          
          for (t in 1:G){
            weighted_alpha[[t]] <- weighted_alpha[[t]] * apply(safelog(phiD_list[[t]]) + (weighted_alpha[[t]] %*% diag(coef[1:K+K])), 2, softmax)
          }

          alpha_agg <- Reduce('rbind', lapply(weighted_alpha, function(x) diag(colSums(x))))
          ###### LASSO PORTION ######
          mod <- NULL
          inner_count <- 1
          while(is.null(mod)) {
            nlambda <- 250
            mod <- tryCatch(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                           offset=offset_agg, standardize=FALSE,
                                           intercept=FALSE,
                                           lambda.min.ratio=0.001,
                                           nlambda=nlambda, alpha=1,lower.limits = low, upper.limits = upp,
                                           maxit=10000, thresh=1e-05, 
                                           penalty.factor = penalty.factor),
                            warning=function(w) return(NULL),
                            error=function(e) stop(e))
            #if it didn't converge, increase nlambda paths by 50%
            if(is.null(mod)) nlambda <- nlambda + floor(50*nlambda)
            inner_count <-inner_count + 1
            #if it didn't converge after 10, just take the estimates as is
            if (inner_count == 10) {
              mod <- suppressWarnings(glmnet::glmnet(x=cbind(x1_agg, alpha_agg), y=y_agg, family="poisson",
                                                     offset=offset_agg, standardize=FALSE,
                                                     intercept=FALSE,
                                                     lambda.min.ratio=0.001,
                                                     nlambda=nlambda, alpha=1,lower.limits = low, upper.limits = upp,
                                                     maxit=10000, thresh=1e-05, 
                                                     penalty.factor = penalty.factor))
            }
          }
          dev <- (1-mod$dev.ratio)*mod$nulldev
          ic <- dev + 2*mod$df ## AIC
          penalty <- which.min(ic)
          coef <- mod$beta[,penalty] #return coefficients
          ###### END LASSO PORTION ######
          
                    
          delta <- coef - coef_old 
          delta_converge <- mean(abs((delta))/abs(0.1+coef_old))
          
          if(delta_converge <= 1e-4) {break}
          kk <- kk + 1
          
        }
        
        
      } 
      
      list(coef = coef)
      
    }

    phiD_list <- lapply(split(phiD, group), matrix, ncol = K) 
    alpha_s <- lapply(split(alphaS, group), matrix, ncol = K) ## numGroups by topics
    G <- length(alpha_s)
    
    # Register the parallel backend
    # cl <- makeCluster(detectCores() - 1)
    # registerDoParallel(cl)
    
    # Use foreach for parallel processing
    i = 1
    results <- foreach(i = 1:V, .export = c('softmax', 'safelog', 'lse', 'alpha_agg')) %dopar% {
      process_word(i)
    }
    # Stop the cluster when done
    # stopCluster(cl)

    # results <- mclapply(1:V, process_word, mc.cores = detectCores() - 1)
    
    # cl <- makeCluster(detectCores() - 1)
    
    # Export necessary variables and functions to the cluster
    # results <- parLapply(cl, 1:V, process_word)
    # stopCluster(cl)
    
    for (v in 1:V) {
      # kappa$kappa_t[v,] <- results[[v]]$coef[1:K]
      # kappa$kappa_s[v,] <- results[[v]]$coef[1:K+K]
      kappa$kappa_t[v,] <- (kappa$kappa_t[v,] + results[[v]]$coef[1:K])/2
      kappa$kappa_s[v,] <- (kappa$kappa_s[v,] + results[[v]]$coef[1:K+K])/2
    }
  }
  return(kappa)
}

