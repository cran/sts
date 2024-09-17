#' Heldout Log-Likelihood
#' 
#' Compute the heldout log-likelihood of the STS model
#' 
#' @param mv the baseline log-transformed occurrence rate of each word in the 
#' corpus
#' @param kappa the estimated kappa coefficients
#' @param alpha the estimated alpha values for the corpus
#' @param missing list of which words and documents are in the heldout set
#' 
#' @return expected.heldout is the average of the held-out log-likelihood values 
#' for each document.
#' 
#' @examples 
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' X <- model.matrix(~1+out$meta$treatment + out$meta$pid_rep + 
#' out$meta$treatment * out$meta$pid_rep)[,-1]
#' X_seed <- as.matrix(out$meta$treatment)
#' out <- make.heldout(out$documents, out$vocab)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(X, X_seed, out, numTopics = 3, verbose = FALSE, 
#' parallelize = FALSE, maxIter = 3, initialization = 'anchor')
#' sm <- sample(x=1:length(out$missing$index), 
#' size = length(out$missing$index)*0.8, replace = TRUE)
#' d.h <- list(index = out$missing$index[sm], docs = out$missing$docs[sm])
#' heldoutLikelihood(mv=sts_estimate$mv, kappa=sts_estimate$kappa, 
#' alpha=sts_estimate$alpha, missing=d.h)$expected.heldout
#' @export
heldoutLikelihood <- function (mv, kappa, alpha, missing) 
{
  K <- (1 + ncol(alpha))/2
  
  heldout <- vector(length = length(missing$index))
  ntokens <- vector(length = length(missing$index))
  # beta <- lapply(model$beta$logbeta, exp)
  # bindex <- model$settings$covariates$betaindex[missing$index]
  
  expeta <- exp(cbind(alpha[,1:(K-1)],0))
  theta <- expeta/rowSums(expeta)
  alpha_s <- alpha[,1:K+K-1]
  
  for (i in 1:length(missing$index)) {
    docid <- missing$index[i]
    words <- missing$docs[[i]][1, ]

    full_kappa <- exp(mv + kappa$kappa_t + kappa$kappa_s %*% diag(alpha_s[docid,]))
    beta <- t(apply(full_kappa, 1, function(m) m / colSums(full_kappa)))

    # probs <- model$theta[docid, ] %*% beta[,words]
    probs <- theta[docid, ] %*% t(beta)[,words]
    probs <- rep(probs, missing$docs[[i]][2, ])
    heldout[i] <- mean(log(probs))
    ntokens[i] <- sum(missing$docs[[i]][2,])
  }
  out <- list(expected.heldout = mean(heldout, na.rm = TRUE), 
              doc.heldout = heldout, index = missing$index, ntokens = ntokens)
  return(out)
}
