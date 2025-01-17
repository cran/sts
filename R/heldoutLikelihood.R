#' Compute Heldout Log-Likelihood
#' 
#' Compute the heldout log-likelihood of the STS model
#' 
#' @param object an sts object, typically after applying \code{\link[stm]{make.heldout}}
#' @param missing list of which words and documents are in the heldout set
#' 
#' @return expected.heldout is the average of the held-out log-likelihood values 
#' for each document.
#' 
#' @examples \donttest{
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' out_ho <- make.heldout(out$documents, out$vocab)
#' out_ho$meta <- out$meta
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out_ho, K = 3, maxIter = 2, verbose = FALSE)
#' heldoutLikelihood(sts_estimate, out_ho$missing)$expected.heldout
#' }
#' @export
heldoutLikelihood <- function (object, missing) 
{
  mv <- object$mv
  kappa <- object$kappa
  alpha <- object$alpha
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
