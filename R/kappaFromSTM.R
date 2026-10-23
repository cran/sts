## Initialize the STS word-choice parameters (kappa_t, kappa_s) from an STM
## fitted with initializationVar as its content covariate
## (initialization = "stm").
##
## In the STM with content covariate levels a = 1, ..., A and topic-covariate
## interactions, the log-rate of word v in topic k for documents in level a is
##   m_v + kappa_topic[v, k] + kappa_level[v, a] + kappa_int[v, k, a],
## where m is the same baseline log word frequency as sts's mv. In the STS
## model it is mv_v + kappa_t[v, k] + kappa_s[v, k] * alpha_s, and alpha_s is
## initialized to the centered initializationVar, x_a - mean(x).
##
## For each topic, kappa_s is the slope and (kappa_t - kappa_topic) the
## intercept of a least-squares fit of the STM content effect
## (kappa_level + kappa_int) on the centered level values, with each level
## weighted by its number of documents. With two levels (e.g., a binary
## initializationVar) this reproduces the STM's two topic-word distributions
## exactly; with more levels it is their best linear approximation.
##
## The STM stores the content-model coefficients in beta$kappa$params as
## K topic effects, then A level effects, then K*A interactions ordered with
## topic varying fastest: params[[K + A + (a - 1) * K + k]].
kappaFromSTM <- function(stmfit, K, V, vals, group) {
  params <- stmfit$beta$kappa$params
  A <- (length(params) - K) / (K + 1)
  numGroups <- length(vals)
  if (is.null(params) || A < 2 || A != numGroups) {
    stop("initialization = \"stm\" requires an STM fit with initializationVar as its ",
         "content covariate (with topic-covariate interactions). ",
         sprintf("The STM has %s content covariate levels but initializationVar has %i distinct values.",
                 if (is.null(params) || A < 1 || A != round(A)) "no" else A, numGroups))
  }
  if (any(lengths(params) != V)) {
    stop("The STM vocabulary does not match the corpus vocabulary.")
  }
  w <- tabulate(group, nbins = numGroups)
  w <- w / sum(w)                          # share of documents in each level
  xc <- vals - sum(w * vals)               # centered level values
  sxx <- sum(w * xc^2)
  kappa_t <- matrix(0, V, K)
  kappa_s <- matrix(0, V, K)
  for (k in 1:K) {
    content_k <- vapply(1:A, function(a) params[[K + a]] + params[[K + A + (a - 1) * K + k]],
                        numeric(V))        # V x A matrix of content effects
    content_k <- matrix(content_k, nrow = V)
    kappa_t[, k] <- params[[k]] + as.numeric(content_k %*% w)
    kappa_s[, k] <- as.numeric(content_k %*% (w * xc)) / sxx
  }
  list(kappa_t = kappa_t, kappa_s = kappa_s)
}
