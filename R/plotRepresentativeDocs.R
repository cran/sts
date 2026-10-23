#' Plot Documents that Load Heavily on a Topic
#' 
#' Produces a plot of the text of documents that load most heavily on topics for an STS object 
#' 
#' @param object Model output from sts.
#' @param text.cex Size of the text; Defaults to 1
#' @param width Size of the plotting window; Defaults to 100
#' @examples
#' \donttest{
#' #Examples with the Gadarian Data
#' library("sts")
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

## Plot a set of quotes (text) stacked vertically, separated by dashed lines.
## Adapted from plotQuote() in the stm package (version 1.3.8; MIT license,
## Copyright (c) 2017 Margaret Roberts, Brandon Stewart and Dustin Tingley);
## see inst/COPYRIGHTS.
plotQuote <- function(sentences, width = 30, text.cex = 1, maxwidth = NULL,
                      main = NULL, xlab = "", ylab = "", xlim = NULL,
                      ylim = NULL, ...) {
  xaxt <- "n"
  yaxt <- "n"
  numlines <- c()
  out <- list()
  for (j in 1:length(sentences)) {
    sentence <- sentences[j]
    if (!is.null(maxwidth)) {
      sentence <- strwrap(sentence, width = maxwidth)[1]
    }
    out[[j]] <- paste(strwrap(sentence, width = width), collapse = "\n")
    numlines[j] <- length(strsplit(out[[j]], "\n")[[1]])
  }
  if (is.null(xlim)) xlim <- c(0, 5)
  if (is.null(ylim)) ylim <- c(0, .5 * sum(numlines))
  graphics::plot(c(0, 0), type = "n", xlim = xlim, ylim = ylim, xaxt = xaxt,
                 yaxt = yaxt, xlab = xlab, ylab = ylab, main = main, ...)
  numlines <- c(0, rev(numlines))
  out <- rev(out)
  start <- 0
  for (j in 2:(length(out) + 1)) {
    graphics::text(2.5, start + numlines[j - 1] / 4 + numlines[j] / 4,
                   out[[j - 1]], cex = text.cex)
    start <- start + numlines[j - 1] / 4 + numlines[j] / 4
    if (j != 1 & j != (length(out) + 1)) {
      graphics::lines(c(0, 5), c(sum(numlines[1:j - 1]) / 2 + numlines[j] / 2,
                                 sum(numlines[1:j - 1]) / 2 + numlines[j] / 2),
                      lty = 2)
    }
  }
  invisible(NULL)
}
