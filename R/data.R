#' Gadarian and Albertson data
#'
#' This data set contains variables from Gadarian and Albertson (2014). The
#' experiment had those in the treatment condition write about what made them
#' anxious about immigration. The control condition just had subjects write
#' about immigration. The data set is also distributed with the stm package
#' (Roberts, Stewart and Tingley; MIT license) and is included here so that
#' the examples do not require stm.
#'
#' @name gadarian
#' @docType data
#' @usage gadarian
#' @format A data frame with 341 observations on the following 4 variables.
#' \describe{
#' \item{\code{MetaID}}{A numeric vector containing identification
#' numbers; not used for analysis}
#' \item{\code{treatment}}{A numeric vector indicating treatment condition}
#' \item{\code{pid_rep}}{A numeric vector of party identification}
#' \item{\code{open.ended.response}}{A character vector of the subject's open
#' ended response}
#' }
#' @source Gadarian, Shana Kushner, and Bethany Albertson. "Anxiety,
#' immigration, and the search for information." Political Psychology 35.2
#' (2014): 133-164.
#'
#' Roberts, Margaret E., Brandon M. Stewart, Dustin Tingley, Christopher Lucas,
#' Jetson Leder-Luis, Shana Kushner Gadarian, Bethany Albertson, and David G.
#' Rand.  "Structural Topic Models for Open-Ended Survey Responses." American
#' Journal of Political Science 58, no 4 (2014): 1064-1082.
#' @keywords datasets
#' @examples
#' head(gadarian)
#' temp <- textProcessor(documents = gadarian$open.ended.response,
#'                       metadata = gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
"gadarian"
