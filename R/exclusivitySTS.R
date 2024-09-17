#' Exclusivity
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
#' @param beta the beta probability  matrix (topic-word distributions) for a given document or alpha-level
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
#' X <- model.matrix(~1+out$meta$treatment + out$meta$pid_rep + 
#' out$meta$treatment * out$meta$pid_rep)[,-1]
#' X_seed <- as.matrix(out$meta$treatment)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(X, X_seed, out, numTopics = 3, verbose = FALSE, 
#' parallelize = FALSE, maxIter = 3, initialization = 'anchor')
#' full_beta_distn <- exp(sts_estimate$mv + sts_estimate$kappa$kappa_t + 
#' sts_estimate$kappa$kappa_s %*% diag(apply(sts_estimate$alpha[,3:5], 2, mean)))
#' full_beta_distn <- t(apply(full_beta_distn, 1, 
#' function(m) m / colSums(full_beta_distn)))
#' exclusivitySTS(full_beta_distn)
#' }
#' @export
exclusivitySTS = function (beta, M = 10, frexw = 0.7) 
{
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
