

whatiscorrect <- function(votematrix) {
  correct.class.vote <- numeric(nrow(votematrix))
  correct.prediction <- logical(nrow(votematrix))
  best.vote          <- character(nrow(votematrix))
  for(i in 1:nrow(votematrix)) {
    correct.class.vote[i] <- votematrix[i, rownames(votematrix)[i]==colnames(votematrix) ]
    correct.prediction[i] <- rownames(votematrix)[i] %in% colnames(votematrix)[votematrix[i,]==max(votematrix[i,])]
   best.vote[i] <- ifelse(correct.prediction[i],
                          rownames(votematrix)[i],
                          colnames(votematrix)[which.max(votematrix[i,])])
  }
  return(list(correct.class.vote=correct.class.vote,
              correct.prediction=correct.prediction,
              best.vote=best.vote))
}
