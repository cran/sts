#' Print estimated regression tables
#' 
#' Prints estimated regression tables from estimateRegnTables()
#' 
#' @param x the estimated regression tables from estimateRegnTables()
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
#' X <- model.matrix(~1+out$meta$treatment + out$meta$pid_rep + 
#' out$meta$treatment * out$meta$pid_rep)[,-1]
#' X_seed <- as.matrix(out$meta$treatment)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(X, X_seed, out, numTopics = 3, verbose = FALSE, 
#' parallelize = FALSE, maxIter = 3, initialization = 'anchor')
#' regn_tables <- estimateRegnTables(sts_estimate$alpha, mu = sts_estimate$mu, 
#' sigma = sts_estimate$sigma, X = X)
#' printRegnTables(regn_tables)
#' }
#' @export
printRegnTables <- function(x, digits = max(3L, getOption("digits") - 3L), 
                             signif.stars = getOption("show.signif.stars"), ...) {
#  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
#      "\n\n", sep = "")
  K <- length(x)/2
  for(i in 1:length(x)) {
    if (i <= K) {
      cat(sprintf("\nTopic %d alpha_p:\n", i))
    } else {
      cat(sprintf("\nTopic %d alpha_s:\n", i - K))
    }
    cat("\nCoefficients:\n")
    coefs <- x[[i]]
    stats::printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                        na.print = "NA", ...)
    cat("\n")
  }
  invisible(x)
}
