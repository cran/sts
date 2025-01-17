#' Compute Exclusivity
#' 
#' Calculate an exclusivity metric for an STS model.
#' 
#' Roberts et al 2014 proposed an exclusivity measure to help with topic model 
#' selection. 
#' 
#' The exclusivity measure includes some information on word frequency as well.  
#' It is based on the FREX
#' labeling metric (see Roberts et al. 2014) with the weight set to .7 in 
#' favor of exclusivity by default.
#'  
#' @param object Model output from sts
#' @param M the number of top words to consider per topic
#' @param frexw the frex weight
#' 
#' @return a numeric vector containing exclusivity for each topic
#' 
#' @references 
#' Mimno, D., Wallach, H. M., Talley, E., Leenders, M., and 
#' McCallum, A. (2011, July). "Optimizing semantic coherence in topic models." 
#' In Proceedings of the Conference on Empirical Methods in 
#' Natural Language Processing (pp. 262-272). Association for 
#' Computational Linguistics. Chicago
#' 
#' Bischof and Airoldi (2012) "Summarizing topical content with word frequency 
#' and exclusivity" In Proceedings of the International Conference on 
#' Machine Learning.
#' 
#' Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., 
#' Gadarian, S., Albertson, B., et al. (2014). 
#' "Structural topic models for open ended survey responses." American Journal 
#' of Political Science, 58(4), 1064-1082.
#' @examples  \donttest{
#' #An example using the Gadarian data from the stm package.  
#' # From Raw text to fitted model using textProcessor() which leverages the 
#' # tm Package
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' topicExclusivity(sts_estimate)
#' }
#' @export
topicExclusivity = function (object, M = 10, frexw = 0.7) 
{
  K <- (1 + ncol(object$alpha))/2
  beta <- exp(object$mv + object$kappa$kappa_t + object$kappa$kappa_s %*% diag(apply(object$alpha[,1:K+K-1], 2, mean)))
  beta <- t(apply(beta, 1, function(m) m / colSums(beta)))
    
  w <- frexw
  tbeta <- beta
  s <- rowSums(tbeta)
  mat <- tbeta/s
  ex <- apply(mat, 2, rank)/nrow(mat)
  fr <- apply(tbeta, 2, rank)/nrow(mat)
  frex <- 1/(w/ex + (1 - w)/fr)
  index <- apply(tbeta, 2, order, decreasing = TRUE)[1:M, ]
  out <- vector(length = ncol(tbeta))
  for (i in 1:ncol(frex)) {
    out[i] <- sum(frex[index[, i], i])
  }
  out
}
