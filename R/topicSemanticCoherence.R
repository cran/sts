#' Compute Semantic Coherence
#' 
#' Calculates semantic coherence for an STS model.
#' 
#' @param object Model output from sts
#' @param corpus The document term matrix to be modeled in a sparse term count matrix with one row
#' per document and one column per term. The object must be a list of with each element 
#' corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is the same format in the \code{\link[stm]{stm}} package. 
#' @param M the number of top words to consider per topic
#' 
#' @return a numeric vector containing semantic coherence for each topic
#' 
#' @examples \donttest{
#' #An example using the Gadarian data from the stm package.  From Raw text to 
#' # fitted model using textProcessor() which leverages the tm Package
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2, verbose = FALSE)
#' topicSemanticCoherence(sts_estimate, out)
#' }
#' @export
topicSemanticCoherence = function (object, corpus, M = 10) 
{
  semCoh1beta <- function(mat, M, beta){
    #Get the Top N Words
    top.words <- apply(beta, 1, order, decreasing=TRUE)[1:M,]
    wordlist <- unique(as.vector(top.words))
    mat <- mat[,wordlist]
    mat$v <- ifelse(mat$v>1, 1,mat$v) #binarize
    
    #do the cross product to get co-occurrences
    cross <- slam::tcrossprod_simple_triplet_matrix(t(mat))
    
    #create a list object with the renumbered words (so now it corresponds to the rows in the table)
    temp <- match(as.vector(top.words),wordlist)
    labels <- split(temp, rep(1:nrow(beta), each=M)) ### CHECK
    
    #Note this could be done with recursion in an elegant way, but let's just be simpler about it.
    sem <- function(ml,cross) {
      m <- ml[1]; l <- ml[2]
     # log(.01 + cross[m,l]) - log(cross[l,l] + .01)  ## THIS LINE MATCHES THE STM PACKAGE
      log(1 + cross[m,l]) - log(cross[l,l]) ## THIS LINE WILL MATCH THE ORIGINAL PAPER BY MIMNO ET AL (Optimizing Semantic Coherence in Topic Models)
    }
    result <- vector(length=nrow(beta))
    for(k in 1:nrow(beta)) {
      grid <- expand.grid(labels[[k]],labels[[k]])
      colnames(grid) <- c("m", "l") #corresponds to original paper
      grid <- grid[grid$m > grid$l,]
      calc <- apply(grid,1,sem,cross)
      result[k] <- sum(calc)
    }
    return(result)
  }
  
  K <- (1 + ncol(object$alpha))/2
  beta <- exp(object$mv + object$kappa$kappa_t + object$kappa$kappa_s %*% diag(apply(object$alpha[,1:K+K-1], 2, mean)))
  beta <- t(apply(beta, 1, function(m) m / colSums(beta)))
  
  documents <- corpus$documents
  vocab <- corpus$vocab
  
  args <- asSTMCorpus(documents)
  documents <- args$documents
  beta <- log(beta)
  # triplet <- doc.to.ijv(documents)
  # mat <- slam::simple_triplet_matrix(triplet$i, triplet$j, 
  #                                    triplet$v, ncol = model$settings$dim$V)
  mat <- convertCorpus(documents, vocab, type="slam")
  result = semCoh1beta(mat, M, beta = t(beta))
  return(result)
}
