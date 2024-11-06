#' Function for printing top words that load heavily on each topic
#' 
#' Prints the top words for each document for low, average, and high levels of sentiment-discourse
#' 
#' @param object Model output from sts
#' @param n number of words to print to console for each topic
#' @param lowerPercentile Percentile to calculate a representative negative sentiment document.
#' @param upperPercentile Percentile to calculate a representative positive sentiment document.
#' @examples
#' \donttest{
#' #Examples with the Gadarian Data
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' printTopWords(sts_estimate)
#' }
#' @export printTopWords
#' @export
printTopWords = function(object, n = 10, lowerPercentile = 0.05, upperPercentile = 0.95) {
  mv <- object$mv
  kappa.est <- object$kappa
  alpha.est <- object$alpha
  K <- (ncol(object$sigma)+1)/2
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, mean))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:n]
    object$vocab[windex]
  })
  labs <- apply(topwords_topic, 2, function(x) paste(x,collapse=", "))
  toprint <- sprintf("Topic %i Avg sentiment-discourse: %s \n", 1:length(labs), labs)
  cat(toprint)		
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, quantile, upperPercentile))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:n]
    object$vocab[windex]
  })
  if (length(topwords_topic) > 0) {
    if (is.list(topwords_topic)) {
      labs <- lapply(topwords_topic, function(x) paste(x,collapse=", "))            
    } else{
      labs <- apply(topwords_topic, 2, function(x) paste(x,collapse=", "))
    }
  } else {
    labs <- rep("no positive words", K)
  }
  toprint <- sprintf("Topic %i Positive sentiment-discourse: %s \n", 1:length(labs), labs)
  cat(toprint)
  
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, quantile, lowerPercentile))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:n]
    object$vocab[windex]
  })
  if (length(topwords_topic) > 0) {
    if (is.list(topwords_topic)) {
      labs <- lapply(topwords_topic, function(x) paste(x,collapse=", "))            
    } else{
      labs <- apply(topwords_topic, 2, function(x) paste(x,collapse=", "))
    }
  } else {
    labs <- rep("no negative words", K)
  }
  toprint <- sprintf("Topic %i Negative sentiment-discourse: %s \n", 1:length(labs), labs)
  cat(toprint)
  
  # return(1)        
}
