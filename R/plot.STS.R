#' Function for plotting STS objects
#' 
#' Produces a plot of the most likely words and their probabilities for each topic for different levels 
#' of sentiment for an STS object.  
#' 
#' @param x Model output from sts.
#' @param n Sets the number of words used to label each topic.  In perspective
#' plots it approximately sets the total number of words in the plot. 
#' n must be greater than or equal to 2
#' @param topics Vector of topics to display. Defaults to all topics.
#' @param lowerPercentile Percentile to calculate a representative negative sentiment document.
#' @param upperPercentile Percentile to calculate a representative positive sentiment document.
#' @param \dots Additional parameters passed to plotting functions.

#' @examples
#' \donttest{
#' #Examples with the Gadarian Data
#' library("tm"); library("stm"); library("sts")
#' temp<-textProcessor(documents=gadarian$open.ended.response,
#' metadata=gadarian, verbose = FALSE)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose = FALSE)
#' out$meta$noTreatment <- ifelse(out$meta$treatment == 1, -1, 1)
#' ## low max iteration number just for testing
#' sts_estimate <- sts(~ treatment*pid_rep, ~ noTreatment, out, K = 3, maxIter = 2)
#' plot(sts_estimate)
#' plot(sts_estimate, n = 10, topic = c(1,2))
#' }
#' @export
plot.STS <- function(x, 
                     n=10, topics=NULL,
                     lowerPercentile = 0.05, upperPercentile = 0.95,
                     ...){
  numWords <- n

  alpha.est <- x$alpha
  kappa.est <- x$kappa
  mu <- x$mu
  sigma <- x$sigma
  mv <- x$mv
  K <- (1+ncol(alpha.est))/2
  
    
  # full_beta_distn is calculated often in this script -- it is the beta probability distribution in eqn 2 of the paper
  full_beta_distn <- exp(mv + kappa.est$kappa_t + kappa.est$kappa_s  %*% diag( apply(alpha.est[,1:K+(K-1)], 2, mean)))  
  full_beta_distn <- t(apply(full_beta_distn, 1, function(m) m / colSums(full_beta_distn)))
  df <- data.frame(index = numWords:1, topic = 1, words = x$vocab[order(full_beta_distn[,1], decreasing = TRUE)[1:numWords]], 
                   kappa = sort(full_beta_distn[,1], decreasing = TRUE)[1:numWords])
  for (k in 2:K) {
    df <- rbind(df, data.frame(index = numWords:1, topic = k, words = x$vocab[order(full_beta_distn[,k], decreasing = TRUE)[1:numWords]],
                               kappa = sort(full_beta_distn[,k], decreasing = TRUE)[1:numWords])
    )
  }
  df$Model <- sprintf("STS: Average")
  
  full_beta_distn <- exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+(K-1)], 2, quantile, lowerPercentile)))
  full_beta_distn <- t(apply(full_beta_distn, 1, function(m) m / colSums(full_beta_distn)))
  df2 <- data.frame(index = numWords:1, topic = 1, words = x$vocab[order(full_beta_distn[,1], decreasing = TRUE)[1:numWords]], 
                    kappa = sort(full_beta_distn[,1], decreasing = TRUE)[1:numWords])
  for (k in 2:K) {
    df2 <- rbind(df2, data.frame(index = numWords:1, topic = k, words = x$vocab[order(full_beta_distn[,k], decreasing = TRUE)[1:numWords]],
                                 kappa = sort(full_beta_distn[,k], decreasing = TRUE)[1:numWords])
    )
  }
  df2$Model <- "STS: Negative"
  
  full_beta_distn <- exp(mv + kappa.est$kappa_t + kappa.est$kappa_s %*% diag(apply(alpha.est[,1:K+(K-1)], 2, quantile, upperPercentile)))
  full_beta_distn <- t(apply(full_beta_distn, 1, function(m) m / colSums(full_beta_distn)))
  df3 <- data.frame(index = numWords:1, topic = 1, words = x$vocab[order(full_beta_distn[,1], decreasing = TRUE)[1:numWords]], 
                    kappa = sort(full_beta_distn[,1], decreasing = TRUE)[1:numWords])
  for (k in 2:K) {
    df3 <- rbind(df3, data.frame(index = numWords:1, topic = k, words = x$vocab[order(full_beta_distn[,k], decreasing = TRUE)[1:numWords]],
                                 kappa = sort(full_beta_distn[,k], decreasing = TRUE)[1:numWords])
    )
  }
  df3$Model <- "STS: Positive"
  
  df <- rbind(df, df2, df3)
  
  
  if (!is.null(topics)) {
    df <- df[df$topic %in% topics,]
  }
  
  df$topic <- paste0("Topic ", df$topic)
  df$cols <- sign(df$kappa)
  # df$index2 <- df$index
  index <- df$index
  words <- df$words

  df$topic <- factor(df$topic, levels = paste0("Topic ", 1:10))
  
  df$Model <- factor(df$Model, levels = c("STS: Average", "STS: Negative", "STS: Positive"))
  
  final_plot <- ggplot(df, aes(index, kappa, label = words)) + geom_vline(xintercept = 0, colour = "grey") + 
          geom_col(fill = "lightblue") + 
          geom_text(check_overlap = FALSE, size = 3, hjust = "inward", angle = 0, colour = "darkred") + 
          theme_bw() + xlab("") + ylab("Probability") + 
          facet_grid(Model ~ topic) +  theme(axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),                                    
                                             axis.ticks.y=element_blank(),
                                             axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + 
          guides(colour='none')+ coord_flip()
  
  # print(final_plot)
  return(final_plot)

}