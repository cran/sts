#' Estimate approximate ELBO 
#' 
#' Estimates the approximate ELBO needed for the Variational E-step in C++
#' 
#' @param alpha_d the estimated alpha variables for the given document
#' @param Sigma_Inv the inverse covariance matrix 
#' @param kappa_t the estimated kappa_t coefficients
#' @param kappa_s the estimated kappa_s coefficients
#' @param mu_d the fitted values based on document-level variables * estimated Gamma for each document.
#' @param doc the sparse matrix representation of the document, with two rows, and columns equal to the number of unique
#' vocabulary words in the document.
#' @param V the size of the vocabulary
#' @param mv the baseline log-transformed occurrence rate of each word in the corpus
#' 
#' @return The Hessian matrix
#' @keywords internal
#' @examples \donttest{
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' # for document #1: 
#' aelbo <- lpbdcpp(alpha_d = sts_estimate$alpha[1,], kappa_t=sts_estimate$kappa$kappa_t, 
#' kappa_s=sts_estimate$kappa$kappa_s, Sigma_Inv = sts_estimate$sigma_inv, doc = out$documents[[1]], 
#' V=length(sts_estimate$vocab), mv = sts_estimate$mv, mu_d = sts_estimate$mu[,1])
#' }
#' @export
lpbdcpp <- function(alpha_d, Sigma_Inv, kappa_t, kappa_s, mu_d, doc, V, mv) {
  .Call('_sts_lpbdcpp', alpha_d, Sigma_Inv, kappa_t, kappa_s, mu_d, doc, V, mv)
}
