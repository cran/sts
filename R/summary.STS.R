#' Summary Function for the STS objects
#' 
#' Function to report on the contents of STS objects
#' 
#' Summary prints a short statement about the model and then runs
#' \code{\link{printTopWords}}.
#' 
#' @aliases summary.STS print.STS
#' @param object An STS object.
#' @param \dots Additional arguments affecting the summary
#' @method summary STS
#' @export
summary.STS <- function(object,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     (ncol(object$alpha)+1)/2, 
                     nrow(object$alpha), 
                     nrow(object$kappa$kappa_s))
  cat(toprint)
  cat("Most likely words:\n")
  printTopWords(object)
}

#' @method print STS
#' @export
print.STS <- function(x,...) {
  toprint <- sprintf("A topic model with %i topics, %i documents and a %i word dictionary.\n", 
                     (ncol(x$alpha)+1)/2, 
                     nrow(x$alpha), 
                     nrow(x$kappa$kappa_s))
  cat(toprint)
}