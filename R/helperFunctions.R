##
#A series of fast softmax functions mostly wrappers around matrixStats package functions -- code from the STM package
##
logsoftmax <- function(x) {
  x - lse(x)
}

lse <- function(x) {
 matrixStats::logSumExp(x)
}

row.lse <- function(mat) {
  matrixStats::rowLogSumExps(mat)
}
col.lse <- function(mat) {
  matrixStats::colLogSumExps(mat)
}

softmax <- function(x) {
  exp(x - lse(x))
}

safelog <- function(x, min=-1000) {
  out <- log(x)
  out[which(out< min)] <- min
  out
}

safeexp = function(x, max=1e8) {
  out <- exp(x)
  out[which(out > max)] <- max
  out
}
