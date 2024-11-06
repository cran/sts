#' Function for plotting documents that load heavily on a topic
#' 
#' Produces a plot of the text of documents that load most heavily on topics for an STS object 
#' 
#' @param object Model output from sts.
#' @param text.cex Size of the text; Defaults to 1
#' @param width Size of the plotting window; Defaults to 100
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
#' docs <- findRepresentativeDocs(sts_estimate, out$meta$open.ended.response, topic = 3, n = 1)
#' plotRepresentativeDocs(docs, text.cex = 0.7, width = 100)
#' }
#' @export 
plotRepresentativeDocs = function(object, text.cex=1, width=100) {
  create_text = function(x) {
    if (x == 1) {
      return("bottom-most quartile (0-25%).")
    } else if (x == 2) {
      return("2nd quartile (25-50%).")
    } else if (x == 3) {
      return("3rd quartile (50-75%).")
    } else if (x == 4) {
      return("top-most quartile (75-100%).")
    }
  }
  # titles <- paste0("Estimated Sentiment-Discourse is in the\n",
                   # sapply(object$sentiment_quartile, create_text))
  # plotQuote(object$text, text.cex = text.cex, width = width, main = titles[1])
  plotQuote(object$text, text.cex = text.cex, width = width)
}  
