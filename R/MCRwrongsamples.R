

MCRwrongsamples <- function (x, col.names=names(x),rownames.from.object=TRUE,subgroup=NULL,freq=FALSE)
  {
  if (class(x) != "list") stop("The first argument must be a list.")
  if (!(all(sapply(x,class)=="MCRestimate"))) stop("Each object in your given list must be a member of class 'MCRestimate'.")
  if (!is.null(col.names) & length(x) != length(col.names)) stop("The length of your list and the length of col.names must be equal.")
  correct.prediction <- sapply(x,function(x) x$correct.prediction)
  if (freq) wrong.pred.frq <- sapply(x,function(x) 1- x$correct.class.vote)
  else      wrong.pred.frq <- apply(correct.prediction,1:2,function (x) ifelse(x,"","X"))
  if (is.null(col.names))
    colnames(wrong.pred.frq) <- paste("Method", 1:length(x)) else
    colnames(wrong.pred.frq) <- col.names
  if (rownames.from.object) {ref.name <- x[[1]]$sample.names
                             if (!(all(sapply(x,function(x) all(x$sample.names==ref.name))))) stop("At least one sample name lable does not fit.")
                             rownames(wrong.pred.frq) <- ref.name
                            }
  else rownames(wrong.pred.frq) <- 1:nrow(wrong.pred.frq)
  if (!is.null(subgroup)){ ref.class <- x[[1]]$classes
                           if (!(all(sapply(x,function(x) all(x$classes==ref.class))))) stop("At least one class lable does not fit.")
                           if (! subgroup %in% levels(ref.class)) stop("Your subgroup must be one of the given sample groups.")
                           correct.prediction <- correct.prediction[ref.class==subgroup,,drop=FALSE]
                           wrong.pred.frq     <-wrong.pred.frq[ref.class==subgroup,,drop=FALSE]}
  
  throw.out <- apply(correct.prediction,1,all)
  reduced.correct.prediction <- correct.prediction[!throw.out,,drop=FALSE]
  reduced.wrong.pred.frq     <- wrong.pred.frq[!throw.out,,drop=FALSE]

  the.order <- order(rowSums(reduced.correct.prediction))
  ordered.prediction <- reduced.correct.prediction[the.order,,drop=FALSE]
  the.result.table   <- reduced.wrong.pred.frq[the.order,,drop=FALSE]
  the.result.table[ordered.prediction]  <- ""
  return(the.result.table)}


