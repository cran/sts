#' A Structural Topic and Sentiment-Discourse Model for Text Analysis
#' 
#' This package implements the Structural Topic and Sentiment-Discourse  
#' (STS) model, which allows researchers to estimate topic models with 
#' document-level metadata that determines both topic prevalence and 
#' sentiment-discourse. The sentiment-discourse is modeled as a document-level 
#' latent variable for each topic that modulates the word frequency within a 
#' topic. These latent topic sentiment-discourse variables are controlled by 
#' the document-level metadata. The STS model can be useful for regression 
#' analysis with text data in addition to topic modeling's traditional 
#' use of descriptive analysis.
#' 
#' Function to fit the model: \code{\link{sts}} 
#' 
#' Functions for Post-Estimation: \code{\link{estimateRegns}}
#' \code{\link{topicExclusivity}} \code{\link{topicSemanticCoherence}}
#'  \code{\link{heldoutLikelihood}} \code{\link{plotRepresentativeDocs}}
#'  \code{\link{findRepresentativeDocs}} \code{\link{printTopWords}}
#'   \code{\link{plot.STS}}
#' 
#' 
#' @name sts-package
#' @docType package
#' @author Author: Shawn Mankad and Li Chen
#' 
#' Maintainer: Shawn Mankad <smankad@@ncsu.edu>
#' @seealso \code{\link{sts}}
#' @references Chen L. and Mankad, S. (2024) "A Structural Topic and 
#' Sentiment-Discourse Model for Text Analysis" Management Science.
#' @keywords package
#' 
#' @import glmnet
#' @import matrixStats
#' @import slam
#' @import stm
#' @import Matrix
#' @import mvtnorm
#' @import tm
#' @import ggplot2
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel mclapply
#' @importFrom parallel stopCluster
#' @importFrom stats model.matrix
#' @importFrom stats optim
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @useDynLib sts, .registration = TRUE
NULL