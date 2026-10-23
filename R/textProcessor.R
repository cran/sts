##############################################################################
## Text preparation utilities: textProcessor(), prepDocuments() and
## make.heldout().
##
## These functions reproduce the functions of the same name in the stm package
## so that a corpus can be prepared for sts() without installing stm. They are
## adapted from stm version 1.3.8 (files R/textProcessor.R, R/prepDocuments.R,
## R/heldout.R and R/readCorpus.R), written by Margaret E. Roberts,
## Brandon M. Stewart and Dustin Tingley and released under the MIT license:
##   Copyright (c) 2017 Margaret Roberts, Brandon Stewart and Dustin Tingley
## See inst/COPYRIGHTS for the full notice.
##############################################################################

#' Process a Vector of Raw Texts
#'
#' Function that takes in a vector of raw texts (in a variety of languages) and
#' performs basic operations.  This function is essentially a wrapper around the
#' \pkg{tm} package where various user specified options can be selected. It is
#' a port of the function of the same name in the stm package, so that a corpus
#' can be prepared for \code{\link{sts}} without installing stm.
#'
#' This function is designed to provide a convenient and quick way to process a
#' relatively small volume of texts for analysis with the package. It is
#' designed to quickly ingest data in a simple form like a spreadsheet where
#' each document sits in a single cell. Once the text has been processed by
#' \code{textProcessor}, it is usually passed to \code{\link{prepDocuments}}.
#'
#' The processor always strips extra white space and removes non-printing
#' characters. Other operations (lower-casing, stop word, number and
#' punctuation removal, and stemming) are optional. Stemming uses the Snowball
#' stemmer from the \pkg{SnowballC} package via \code{tm::stemDocument}.
#'
#' @param documents The documents to be processed.  A character vector where
#' each entry is the full text of a document (if passed as a different type it
#' will attempt to convert to a character vector).
#' @param metadata Additional data about the documents.  Specifically a
#' \code{data.frame} or \code{matrix} object with number of rows equal to the
#' number of documents and one column per meta-data type. The column names are
#' used to label the metadata.  The metadata do not affect the text processing,
#' but providing the metadata object ensures that if documents are dropped the
#' corresponding metadata rows are dropped as well.
#' @param lowercase Whether all words should be converted to lower case.
#' Defaults to TRUE.
#' @param removestopwords Whether stop words should be removed using the
#' stop word list in \pkg{tm}. Defaults to TRUE.
#' @param removenumbers Whether numbers should be removed. Defaults to TRUE.
#' @param removepunctuation whether punctuation should be removed.  Defaults to
#' TRUE.
#' @param ucp When TRUE passes \code{ucp=TRUE} to \code{tm::removePunctuation}
#' which removes a more general set of punctuation (the Unicode general category
#' P).  Defaults to \code{FALSE}.
#' @param stem Whether or not to stem words. Defaults to TRUE.
#' @param wordLengths From the \pkg{tm} package. An integer vector of length 2.
#' Words shorter than the minimum word length \code{wordLengths[1]} or longer
#' than the maximum word length \code{wordLengths[2]} are discarded. Defaults
#' to \code{c(3, Inf)}, i.e., a minimum word length of 3 characters.
#' @param sparselevel Removes terms where at least sparselevel proportion of
#' the entries are 0. Defaults to 1 which effectively turns the feature off.
#' @param language Language used for processing. Defaults to English. \code{tm}
#' uses the \code{SnowballC} stemmer which as of version 0.5 supports "danish
#' dutch english finnish french german hungarian italian norwegian portuguese
#' romanian russian spanish swedish turkish".  These can be specified as any on
#' of the above strings or by the three-letter ISO-639 codes.  You can also set
#' language to "na" if you want to leave it deliberately unspecified (see
#' documentation in \code{tm}).
#' @param verbose If true prints information as it processes.
#' @param onlycharacter When TRUE, runs a regular expression substitution to
#' replace all non-alphanumeric characters. These characters can crash
#' textProcessor for some operating systems.  May remove foreign characters
#' depending on encoding. Defaults to FALSE.
#' @param striphtml When TRUE, runs a regular expression substitution to strip
#' html contained within <>.  Defaults to FALSE.
#' @param customstopwords A character vector containing words to be removed.
#' Defaults to NULL which does not remove any additional words.  This function
#' is primarily for easy removal of application specific stopwords.  Note that
#' as with other stopwords these are removed after converting everything to
#' lower case but before removing numbers, punctuation or stemming.  Thus words
#' to be removed should be all lower case but otherwise complete.
#' @param custompunctuation A character vector containing any characters to be
#' removed immediately after standard punctuation removal. This function exists
#' primarily for easy removal of application specific punctuation not caught by
#' the punctuation filter (although in principle it could be used to remove
#' any term). If you pass a single character string starting with an open
#' bracket, it is treated as a regular expression; otherwise the characters are
#' combined into a bracketed character class.
#' @param v1 A logical which defaults to \code{FALSE}.  If set to \code{TRUE}
#' it will use the ordering of operations used in Version 1.0 of the stm
#' package to remove punctuation.  Primarily for replication purposes.
#'
#' @return A list containing
#' \item{documents}{A list of documents in the format used by \code{\link{sts}}.}
#' \item{vocab}{Character vector of vocabulary.}
#' \item{meta}{Data frame or matrix containing the user-supplied metadata for
#' the retained documents.}
#' \item{docs.removed}{Indices of the original documents that were dropped
#' because they contained no words after processing.}
#'
#' @seealso \code{\link{prepDocuments}}
#' @references
#' Ingo Feinerer and Kurt Hornik (2013). tm: Text Mining Package. R
#' package version 0.5-9.1.
#'
#' Ingo Feinerer, Kurt Hornik, and David Meyer (2008). Text Mining
#' Infrastructure in R. \emph{Journal of Statistical Software} 25(5): 1-54.
#'
#' Roberts, M. E., Stewart, B. M., and Tingley, D. (2019). "stm: An R Package
#' for Structural Topic Models." \emph{Journal of Statistical Software}
#' 91(2): 1-40.
#' @examples
#' head(gadarian)
#' temp <- textProcessor(documents = gadarian$open.ended.response,
#'                       metadata = gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' @export
textProcessor <- function(documents, metadata = NULL,
                          lowercase = TRUE, removestopwords = TRUE,
                          removenumbers = TRUE, removepunctuation = TRUE,
                          ucp = FALSE, stem = TRUE, wordLengths = c(3, Inf),
                          sparselevel = 1, language = "en", verbose = TRUE,
                          onlycharacter = FALSE, striphtml = FALSE,
                          customstopwords = NULL, custompunctuation = NULL,
                          v1 = FALSE) {

  # as stm does (through stringr), the text is handled as UTF-8 and bytes
  # that are not valid UTF-8 are replaced by U+FFFD
  documents <- asUTF8(documents)

  if (striphtml) {
    documents <- gsub('<.+?>', ' ', documents)
  }

  # remove non-visible characters
  documents <- gsub("(*UCP)[^[:graph:]]", " ", documents, perl = TRUE)

  if (onlycharacter) {
    documents <- gsub("[^[:alnum:]///' ]", " ", documents)
  }

  if (verbose) cat("Building corpus... \n")
  txt <- tm::VCorpus(tm::VectorSource(documents),
                     readerControl = list(language = language))
  # Apply filters
  txt <- tm::tm_map(txt, tm::stripWhitespace)

  if (lowercase) {
    if (verbose) cat("Converting to Lower Case... \n")
    txt <- tm::tm_map(txt, tm::content_transformer(tolower))
  }

  removePunct <- function(txt) {
    if (removepunctuation) {
      if (verbose) cat("Removing punctuation... \n")
      txt <- tm::tm_map(txt, tm::removePunctuation,
                        preserve_intra_word_dashes = TRUE, ucp = ucp)
    }
    if (!is.null(custompunctuation)) {
      if (verbose) cat("Removing custom punctuation... \n")
      if (length(custompunctuation) == 1 &&
          substr(custompunctuation, 0, 1) == "[") {
        # a single entry starting with an open bracket is used as a regex
        punct_pattern <- custompunctuation
      } else {
        punct_pattern <- sprintf("[%s]", paste0(custompunctuation, collapse = ""))
      }
      txt <- tm::tm_map(txt, tm::content_transformer(function(x, pattern) gsub(pattern, "", x)),
                        punct_pattern)
    }
    txt
  }

  if (!v1) txt <- removePunct(txt)
  if (removestopwords) {
    if (verbose) cat("Removing stopwords... \n")
    txt <- tm::tm_map(txt, tm::removeWords, tm::stopwords(language))
  }
  if (!is.null(customstopwords)) {
    if (verbose) cat("Remove Custom Stopwords...\n")
    txt <- tm::tm_map(txt, tm::removeWords, customstopwords)
  }
  if (removenumbers) {
    if (verbose) cat("Removing numbers... \n")
    txt <- tm::tm_map(txt, tm::removeNumbers)
  }
  # the ordering of operations used in version 1.0 of stm
  if (v1) txt <- removePunct(txt)

  if (stem) {
    if (verbose) cat("Stemming... \n")
    txt <- tm::tm_map(txt, tm::stemDocument, language = language)
  }

  # Make a matrix
  if (verbose) cat("Creating Output... \n")
  dtm <- tm::DocumentTermMatrix(txt, control = list(wordLengths = wordLengths,
                                                    tolower = FALSE))
  if (sparselevel != 1) {
    ntokens <- sum(dtm$v)
    V <- ncol(dtm)
    dtm <- tm::removeSparseTerms(dtm, sparselevel)
    if (ncol(dtm) < V & verbose) {
      message <- sprintf("Removed %i of %i terms (%i of %i tokens) due to sparselevel of %s \n",
                         V - ncol(dtm), V, ntokens - sum(dtm$v), ntokens,
                         sparselevel)
      cat(message)
    }
  }
  # If there is metadata we need to remove the rows of dropped documents
  docindex <- sort(unique(dtm$i))
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    if (nrow(metadata) != length(documents)) {
      stop("metadata must have one row per document.")
    }
    rownames(metadata) <- NULL
    metadata <- metadata[docindex, , drop = FALSE]
  }
  out <- read.slam(dtm)

  # Documents with no remaining words are dropped in the conversion above.
  kept <- (1:length(documents) %in% docindex)
  vocab <- as.character(out$vocab)
  out <- list(documents = out$documents, vocab = vocab, meta = metadata,
              docs.removed = which(!kept))
  class(out) <- "textProcessor"
  out
}

#' @method print textProcessor
#' @export
print.textProcessor <- function(x, ...) {
  toprint <- sprintf("A text corpus with %i documents, and an %i word dictionary.\n",
                     length(x$documents), length(x$vocab))
  cat(toprint)
  invisible(x)
}

#' @method summary textProcessor
#' @export
summary.textProcessor <- function(object, ...) {
  toprint <- sprintf("A text corpus with %i documents, and an %i word dictionary. Use str() to inspect object or see documentation \n",
                     length(object$documents), length(object$vocab))
  cat(toprint)
  invisible(object)
}

## Re-encode a character vector as UTF-8, mirroring how stringi (used by stm
## through stringr) reads strings: strings in the native encoding are taken to
## be UTF-8 in UTF-8 and C/POSIX locales and converted from the native
## encoding otherwise; bytes that are not valid UTF-8 are replaced by U+FFFD.
asUTF8 <- function(x) {
  x <- as.character(x)
  unk <- !is.na(x) & Encoding(x) == "unknown"
  li <- l10n_info()
  if (isTRUE(li[["UTF-8"]]) || !(isTRUE(li[["Latin-1"]]) || isTRUE(li[["MBCS"]]))) {
    y <- x[unk]
    Encoding(y) <- "UTF-8"
    x[unk] <- y
  } else {
    x[unk] <- enc2utf8(x[unk])
  }
  x <- enc2utf8(x)
  bad <- !is.na(x) & !validUTF8(x)
  if (any(bad)) x[bad] <- iconv(x[bad], "UTF-8", "UTF-8", sub = "\uFFFD")
  x
}

## Convert a simple triplet matrix (e.g., a tm DocumentTermMatrix) to the
## list-of-documents format
read.slam <- function(corpus) {
  if (!inherits(corpus, "simple_triplet_matrix")) stop("corpus is not a simple triplet matrix")
  if (inherits(corpus, "TermDocumentMatrix")) {
    non_empty_docs <- which(slam::col_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[, non_empty_docs]$j, corpus[, non_empty_docs]$i,
                            corpus[, non_empty_docs]$v)
    names(documents) <- corpus[, non_empty_docs]$dimnames$Docs
  } else {
    non_empty_docs <- which(slam::row_sums(corpus) != 0)
    documents <- ijv.to.doc(corpus[non_empty_docs, ]$i, corpus[non_empty_docs, ]$j,
                            corpus[non_empty_docs, ]$v)
    names(documents) <- corpus[non_empty_docs, ]$dimnames$Docs
  }
  vocab <- corpus$dimnames$Terms
  list(documents = documents, vocab = vocab)
}

#' Prepare Documents for Analysis
#'
#' Performs several corpus manipulations including removing words and
#' renumbering word indices (to correct for zero-indexing and/or unused words
#' in the vocab vector). It is a port of the function of the same name in the
#' stm package, so that a corpus can be prepared for \code{\link{sts}} without
#' installing stm.
#'
#' The default setting \code{lower.thresh=1} means that words which appear in
#' only one document will be dropped.  This is often advantageous as there is
#' little information about these words but the added cost of including them
#' in the model can be quite large.  In many cases it will be helpful to set
#' this threshold considerably higher.  If the vocabulary is in excess of 5000
#' entries inference can slow quite a bit.
#'
#' If words are removed, the function returns a vector of the original indices
#' for the dropped items.  If it removed documents it returns a vector of doc
#' indices removed. Users with accompanying metadata or texts may want to drop
#' those rows from the corresponding objects.
#'
#' If you have a large corpus, consider setting \code{verbose = FALSE}.
#'
#' @param documents List of documents. For more on the format see
#' \code{\link{sts}}.
#' @param vocab Character vector of words in the vocabulary.
#' @param meta Document metadata.
#' @param lower.thresh Words which do not appear in a number of documents
#' greater than lower.thresh will be dropped and both the documents and vocab
#' files will be renumbered accordingly.  If this causes all words within a
#' document to be dropped, a message will print to the screen and it will also
#' return vector of the documents removed so you can update your meta data as
#' well. See details below.
#' @param upper.thresh As with lower.thresh but this provides an upper bound.
#' Words which appear in at least this number of documents will be dropped.
#' Defaults to \code{Inf} which does no filtering.
#' @param subsample If an integer will randomly subsample (without replacement)
#' the given number of documents from the total corpus before any processing.
#' Defaults to \code{NULL} which provides no subsampling.  Note that the output
#' may have fewer than the number of requested documents if additional
#' processing causes some of those documents to be dropped.
#' @param verbose A logical indicating whether or not to print details to the
#' screen.
#' @return A list containing a new documents and vocab object.
#' \item{documents}{The new documents object for use with \code{\link{sts}}}
#' \item{vocab}{The new vocab object} \item{meta}{The new meta data object}
#' \item{docs.removed}{A set of indices corresponding to documents that were
#' removed, if any.} \item{words.removed}{A set of indices corresponding to
#' words that were removed, if any} \item{tokens.removed}{The total number of
#' tokens removed} \item{wordcounts}{A table giving the the number of documents
#' that each word is found in of the original document set, prior to any
#' removal. This can be passed through a histogram for visual inspection}
#' @seealso \code{\link{textProcessor}}
#' @examples
#' temp <- textProcessor(documents = gadarian$open.ended.response,
#'                       metadata = gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' length(out$vocab)
#' @export
prepDocuments <- function(documents, vocab, meta = NULL,
                          lower.thresh = 1, upper.thresh = Inf,
                          subsample = NULL, verbose = TRUE) {
  if ((is.null(documents) | is.null(vocab))) {
    stop("One of your file inputs has no data")
  }
  if (!inherits(documents, "list")) {
    stop("documents must be a list in the format described in ?sts.")
  }

  if (!is.null(subsample)) {
    index <- sample(1:length(documents), subsample)
    documents <- documents[index]
    if (!is.null(meta)) meta <- meta[index, , drop = FALSE]
  }

  # check that there are no 0 length documents
  len <- unlist(lapply(documents, length))
  if (any(len == 0)) {
    stop("Some documents have 0 length.  Please check input.
          See ?prepDocuments() for more info.")
  }

  triplet <- doc.to.ijv(documents) # this also fixes the zero indexing.
  nms <- names(documents)
  documents <- ijv.to.doc(triplet$i, triplet$j, triplet$v)
  names(documents) <- nms
  docs.removed <- c()

  # Detect Missing Terms
  miss.vocab <- NULL
  vocablist <- sort(unique(triplet$j))
  wordcounts <- tabulate(triplet$j)
  if (length(vocablist) > length(vocab)) {
    stop("Your documents object has more unique features than your vocabulary file has entries.")
  }
  if (length(vocablist) < length(vocab)) {
    if (verbose) cat("Detected Missing Terms, renumbering \n")
    miss.vocab <- vocab[-vocablist]
    vocab <- vocab[vocablist]
    new.map <- cbind(vocablist, 1:length(vocablist))
    documents <- lapply(documents, function(d) {
      nm <- names(d)
      d[1, ] <- new.map[match(d[1, ], new.map[, 1]), 2]
      names(d) <- nm
      return(d)
    })
    wordcounts <- wordcounts[vocablist]
  }

  # Remove words appearing in too few / too many documents
  toremove <- which(wordcounts <= lower.thresh | wordcounts >= upper.thresh)
  keepers <- which(wordcounts > lower.thresh & wordcounts < upper.thresh)
  droppedwords <- c(miss.vocab, vocab[toremove])
  if (length(toremove)) {
    if (verbose) cat(sprintf("Removing %i of %i terms (%i of %i tokens) due to frequency \n",
                             length(toremove), length(wordcounts),
                             sum(wordcounts[toremove]), sum(wordcounts)))
    vocab <- vocab[-toremove]
    remap <- 1:length(keepers)
    for (i in 1:length(documents)) {
      doc <- documents[[i]]
      dockeep <- doc[1, ] %in% keepers
      doc <- doc[, dockeep, drop = FALSE]
      doc[1, ] <- remap[match(doc[1, ], keepers)]
      documents[[i]] <- doc
      if (ncol(doc) == 0) docs.removed <- c(docs.removed, i)
    }
    if (length(docs.removed)) {
      if (verbose) cat(sprintf("Removing %i Documents with No Words \n", length(docs.removed)))
      documents <- documents[-docs.removed]
    }
    toprint <- sprintf("Your corpus now has %i documents, %i terms and %i tokens.",
                       length(documents), length(vocab), sum(wordcounts[keepers]))
    if (verbose) cat(toprint)
  }

  if (!is.null(docs.removed) & !is.null(meta)) {
    meta <- meta[-docs.removed, , drop = FALSE]
  }
  # recast everything as an integer
  documents <- lapply(documents, function(x) matrix(as.integer(x), nrow = 2))
  list(documents = documents, vocab = vocab, meta = meta,
       words.removed = droppedwords, docs.removed = docs.removed,
       tokens.removed = sum(wordcounts[toremove]), wordcounts = wordcounts)
}

#' Heldout Likelihood by Document Completion
#'
#' Tools for making and evaluating heldout datasets.  It is a port of the
#' function of the same name in the stm package, so that heldout data can be
#' created for \code{\link{heldoutLikelihood}} without installing stm.
#'
#' These functions are used to create and evaluate heldout likelihood using the
#' document completion method.  The basic idea is to hold out some fraction of
#' the words in a set of documents, train the model and use the document-level
#' latent variables to evaluate the probability of the heldout portion. See the
#' example for the basic workflow.
#'
#' @param documents the documents to be modeled (see \code{\link{sts}} for
#' format).
#' @param vocab the vocabulary item
#' @param N number of docs to be partially held out
#' @param proportion proportion of docs to be held out.
#' @param seed the seed, set for replicability
#' @return A list with elements \code{documents} and \code{vocab} (the
#' training corpus) and \code{missing}, a list with the indices of the
#' partially held-out documents (\code{index}) and their held-out words
#' (\code{docs}).
#' @seealso \code{\link{heldoutLikelihood}}
#' @examples
#' temp <- textProcessor(documents = gadarian$open.ended.response,
#'                       metadata = gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out_ho <- make.heldout(out$documents, out$vocab, seed = 1)
#' length(out_ho$missing$index)
#' @export
make.heldout <- function(documents, vocab, N = floor(.1 * length(documents)),
                         proportion = .5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  index <- sort(sample(1:length(documents), N))
  pie <- proportion
  missing <- vector(mode = "list", length = N)
  ct <- 0
  for (i in index) {
    ct <- ct + 1
    doc <- documents[[i]]
    if (ncol(doc) < 2) next
    doc <- rep(doc[1, ], doc[2, ])
    # how many tokens to sample? The max ensures at least one is sampled
    nsamp <- max(1, floor(pie * length(doc)))
    ho.index <- sample(1:length(doc), nsamp)
    tab <- tabulate(doc[ho.index])
    missing[[ct]] <- rbind(which(tab > 0), tab[tab > 0])
    tab <- tabulate(doc[-ho.index])
    documents[[i]] <- rbind(which(tab > 0), tab[tab > 0])
  }
  missing <- list(index = index, docs = missing)
  # check the vocab
  indices <- sort(unique(unlist(lapply(documents, function(x) x[1, ]))))

  # renumber if some words or held-out documents are missing
  if (length(indices) != length(vocab) | any(unlist(lapply(missing$docs, is.null)))) {
    newind <- rep(0, length(vocab))
    newind[indices] <- 1:length(indices)
    new.map <- cbind(1:length(vocab), newind)
    missing$docs <- lapply(missing$docs, function(d) {
      d[1, ] <- new.map[match(d[1, ], new.map[, 1]), 2]
      return(d[, d[1, ] != 0, drop = FALSE])
    })
    documents <- lapply(documents, function(d) {
      d[1, ] <- new.map[match(d[1, ], new.map[, 1]), 2]
      return(d[, d[1, ] != 0, drop = FALSE])
    })
    lens <- unlist(lapply(missing$docs, length))
    if (any(lens == 0)) {
      missing$docs <- missing$docs[lens != 0]
      missing$index <- missing$index[lens != 0]
    }
    vocab <- vocab[indices]
  }
  heldout <- list(documents = documents, vocab = vocab, missing = missing)
  heldout$documents <- lapply(heldout$documents, function(x) matrix(as.integer(x), nrow(x), ncol(x)))
  heldout$missing$docs <- lapply(heldout$missing$docs, function(x) matrix(as.integer(x), nrow(x), ncol(x)))
  class(heldout) <- "heldout"
  heldout
}
