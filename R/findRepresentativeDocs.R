#' Function for plotting documents that load heavily on a topic
#' 
#' Extracts documents with the highest prevalence for a given topic
#' 
#' @param object Model output from sts
#' @param corpus_text vector of text documents, usually contained in the output of prepDocuments
#' @param topic a single topic number
#' @param n number of documents to extract
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
#' docs <- findRepresentativeDocs(sts_estimate, out$meta$open.ended.response, topic = 3, n = 4)
#' plotRepresentativeDocs(docs, text.cex = 0.7, width = 100)
#' }
#' @export findRepresentativeDocs
#' @export
findRepresentativeDocs = function(object, corpus_text, topic, n = 3) {
  K <- (1 + ncol(object$alpha))/2
  
  expeta <- exp(cbind(object$alpha[,1:(K-1), drop = FALSE],0))
  theta <- expeta/rowSums(expeta)
  
  ind <- order(theta[,topic], decreasing = TRUE)[1:n]
  
  # Calculate quartiles
  quartiles <- quantile(object$alpha[,K + topic - 1], probs = c(0.25, 0.5, 0.75))
  
  # Function to determine the quartile of a value
  get_quartile <- function(value, corpus, quartiles) {
    if (value <= quartiles[1]) {
      return(1)
    } else if (value <= quartiles[2]) {
      return(2)
    } else if (value <= quartiles[3]) {
      return(3)
    } else {
      return(4)
    }
  }
  
  titles <- sapply(object$alpha[ind,K + topic - 1], get_quartile, quartiles = quartiles)
  
  data.frame(text = corpus_text[ind], sentiment_quartile = titles, index = ind)
}  
