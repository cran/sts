# helper function to print the top words to console
print.topWords = function(object) {
  mv <- object$mv
  kappa.est <- object$kappa
  alpha.est <- object$alpha
  K <- (ncol(object$sigma)+1)/2
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, mean))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:10]
    object$vocab[windex]
  })
  labs <- apply(topwords_topic, 2, function(x) paste(x,collapse=", "))
  toprint <- sprintf("Topic %i Avg alpha: %s \n", 1:length(labs), labs)
  cat(toprint)		
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, quantile, 0.9))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:10]
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
  toprint <- sprintf("Topic %i Positive alpha: %s \n", 1:length(labs), labs)
  cat(toprint)
  
  
  topwords_topic <- apply(exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+K-1], 2, quantile, 0.1))), 2, function(x) {
    windex <- order(x,decreasing=TRUE)[1:10]
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
  toprint <- sprintf("Topic %i Negative alpha: %s \n", 1:length(labs), labs)
  cat(toprint)
  
  # return(1)        
}
