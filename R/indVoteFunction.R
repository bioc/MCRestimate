MCRindError <- function(MCRe,perGroup=FALSE){
  if(class(MCRe)!="MCRestimate")
    stop("You have to specify an object of class MCRestimate")
  if(!("indVotes" %in% names(MCRe)))
    stop("Your object does not contain the individual MCR errors. Object was created with an old version of MCRestimate")
  M <- MCRe$"indVotes"
  
  if(perGroup){
    Groups <- colnames(MCRe$"votes")
    WDH    <- dim(M)[3]
    res    <- matrix(0, nrow=2 * (length(Groups)+1), ncol=WDH)
    rownames(res) <- c (Groups, "All", paste(Groups,"fraction",sep="_"), "All_fraction")
    for(i in 1:WDH){
     ma <- M[,,i]
     rownames(ma) <- rownames(MCRe$"votes")
     colnames(ma) <- Groups
     for(j in 1:length(Groups)){
      H    <- Groups[j]
      mi   <- ma[rownames(ma)==H,,drop=FALSE]
      Resu <- whatiscorrect(mi)
      Value <- sum(!Resu$correct.prediction)
      res[H,i] <- Value
      res[paste(H,"fraction",sep="_"),i] <- Value/nrow(mi)
     } 
     Resu  <- whatiscorrect(ma)
     Value <- sum(!Resu$correct.prediction)
     res["All",i]          <- Value
     res["All_fraction",i] <- Value/nrow(ma) 
    }    
  }else{
  res <- c()
  for(i in 1:dim(M)[3]){
    ma <- M[,,i]
    rownames(ma) <- rownames(MCRe$"votes")
    colnames(ma) <- colnames(MCRe$"votes")
    Resu <- whatiscorrect(ma)
    res <- c(res,sum(!Resu$correct.prediction))
  }
  }
  return(res)
}



plotIndGroupVotes <- function(MCRe, PvD= 0.5, dotCol="red", errCol="black", xlab="",ylab="# misclassified samples (mean + SD)",...) { 
  Data   <- MCRindError(MCRe,perGroup=TRUE)
  Groups <- c(colnames(MCRe$"votes"),"All")
  Values <- c()
  Errors <- c()
  MeanFrac <- c()
  for(j in 1:length(Groups)){
   H        <- Groups[j]
   A        <- Data[H,]
   B        <- Data[paste(H,"fraction", sep="_"),]
   MeanFrac <- c(MeanFrac, round(mean(B),2))
   Values   <- c(Values,mean(A))
   Errors   <- c(Errors, sd(A))
 }
  Upper <- Values + Errors
  Lower <- pmax(0,Values - Errors)  
  rangeUp   <- max(Upper) + 1
  rangeDown <- max(0, (min(Lower) - 1))

plot(1:length(Values),Values, xlim=c(0, length(Values) + 1),axes=FALSE, ylim=c(rangeDown, rangeUp),ylab=ylab, xlab=xlab,col=dotCol,...)
axis(2)
box()
for(x in 1:length(Values)){  
  w <- 0.1
  segments(x,Lower[x] , x, Upper[x],col=errCol) #senkrecht
  segments(x - w, Lower[x], x + w, Lower[x],col=errCol) # wagerecht
  segments(x - w, Upper[x], x + w, Upper[x],col=errCol) # wagerecht
  text(x+PvD,Values[x],labels=paste(MeanFrac[x] *100, "%"))
}
points(1:length(Values),Values,col=dotCol,pch=16)
axis(1, at=1:length(Values), labels=Groups,las=2)
}
