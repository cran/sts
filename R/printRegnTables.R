#' Print Estimated Regression Tables
#' 
#' Prints estimated regression tables from estimateRegnTables()
#' 
#' @param x the estimated regression tables from estimateRegnTables()
#' @param topics Vector of topics to display. Defaults to all topics.
#' @param digits minimum number of significant digits to be used for most numbers.
#' @param signif.stars logical; if TRUE, P-values are additionally encoded 
#' visually as ‘significance stars’ in order to help scanning of long 
#' coefficient tables. It defaults to the show.signif.stars slot of options.
#' @param ... other arguments suitable for stats::printCoefmat()
#' 
#' @return Prints estimated regression tables from estimateRegnTables() to console
#' 
#' @examples \donttest{
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' regns <- estimateRegns(sts_estimate, ~treatment*pid_rep, out)
#' printRegnTables(x = regns)
#' }
#' @export
printRegnTables <- function(x, topics = NULL, digits = max(3L, getOption("digits") - 3L), 
                             signif.stars = getOption("show.signif.stars"), ...) {
#  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
#      "\n\n", sep = "")
  K <- length(x)/2
  if (is.null(topics)) {
    prevSentIndex <- 1:length(x)
  } else {
    prevSentIndex <- c(topics, topics + K)
  }
  for(i in prevSentIndex) {
    if (i <= K) {
      cat(sprintf("\nTopic %d prevalence:\n", i))
    } else {
      cat(sprintf("\nTopic %d sentiment:\n", i - K))
    }
    cat("\nCoefficients:\n")
    coefs <- x[[i]]
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                        na.print = "NA", ...)
    cat("\n")
  }
  invisible(x)
}
