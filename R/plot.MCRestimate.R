

plot.MCRestimate <- function(x,
                             class.factor=NULL,
                             rownames.from.object=FALSE,
                             sample.order=TRUE,
                             legend=FALSE,
                             mypalette=NULL,
                             shading=NULL,
                             xlab="Sample ID",
                             ylab="Frequency of correct classification",
                             cex.axis=1,...)
{ require(RColorBrewer)
  if("MCRestimate" %in% class(x))
    { vote.matrix        <- x$votes
      class.factor       <- x$classes
      sample.names       <- x$sample.names
    }
else if(is.matrix(x))
  { vote.matrix      <- x
    sample.names     <- rownames(vote.matrix)
    if (is.null(sample.names))
      sample.names <- 1:nrow(vote.matrix)

    rownames(vote.matrix) <- class.factor
    x                     <- whatiscorrect(vote.matrix)
    
  } else {
    stop("'x' must be a matrix or an object of class MCRestimate")
  }
  correct.prediction <- x$correct.prediction
  correct.class.vote <- x$correct.class.vote
  
  stopifnot(is.factor(class.factor))
  if(!identical(colnames(vote.matrix), levels(class.factor)))
    stop("The column names of the matrix must be the same as 'levels(class.factor)'")
  
  nclasses <- nlevels(class.factor)
  nn      <- length(class.factor)
  
  if (nn!=nrow(vote.matrix))
    stop ("length(class.factor) and the number of rows of the matrix should be the same")
  
  vote.annotation.frame <- data.frame(class.factor,correct.prediction,correct.class.vote, sample.names)
  
  stopifnot(is.logical(sample.order), is.logical(rownames.from.object), is.logical(legend))
  
  if (sample.order)
    {order.index           <- order(class.factor)
     vote.matrix           <- vote.matrix[order.index,]
     vote.annotation.frame <- vote.annotation.frame[order.index,]
   }
  else
    order.index <- 1:nn
  
  ## you get the two highest votes for every sample
  
  two.votes          <- t(apply(vote.matrix,1,sort,decreasing=TRUE)[1:2,])  
  vote.table         <- cbind(two.votes, vote.annotation.frame$correct.class.vote)
  
  ## if a legend is wanted then there is a split of the plot region
  if(legend)
    {mai.save=par("mai")
     a <- mai.save; a[1]<- 0
     par(mai=a)    
     layout(c(1,2),height=c(5,1))
   }
  
  ## The plot is different depending on whether there are two or more classes

  red.color <- brewer.pal(8,"Reds")[6]
  blue.color <- brewer.pal(5,"Blues")[5]
  blue.triagle.color <- brewer.pal(3,"Blues")[2]
  plot(x=0,xlim=c(0,nn),ylim=c(0,1),type="n",axes=FALSE,xlab=xlab,ylab=ylab,...)
  is.correct <- vote.annotation.frame$correct.prediction
  if (nclasses==2) {
    points(1:nn, vote.annotation.frame$correct.class.vote, pch=c(17,19)[1+is.correct],
           col=c(red.color, blue.color)[1+is.correct])
    abline(h=0.5,col="grey")
  } else {
    points(1:nn, two.votes[,1], pch=c(24,19)[1+is.correct], col=c(blue.triagle.color,blue.color)[1+is.correct])
    points(which( is.correct), vote.table[ is.correct, 2], col="gray")
    points(which(!is.correct), vote.table[!is.correct, 3], col=red.color,pch=17)
  }

  
  ## if the argument names is not NULL then the rownames are used to label the x-axis
  if(rownames.from.object) axis(1,at=1:nn,labels=as.character(vote.annotation.frame$sample.names),las=2,cex.axis=cex.axis)
  else                     axis(1,at=1:nn,labels=order.index,cex.axis=cex.axis)
  
  axis(2, las=2,cex.axis=cex.axis)
  
  ## the stripe for the different groups is plotted
  if(is.null(mypalette)) mypalette <- brewer.pal(nclasses+3,"YlGn")[2:(nclasses+1)]
  for (i in 1:nn)
    {color = mypalette[vote.annotation.frame$class.factor[i]]
     if(is.null(shading)){
       dens <- NULL
     }else{
       dens <- shading * (as.numeric(vote.annotation.frame$class.factor[i]))}
     rect(i-0.5,-0.5,i+0.5,-0.015, density=dens,col=color,border=color)
     abline(v=i,lty=2,col="grey")
   }
  
  ## the legend is plotted (if there should be one) 
  if(legend)
    {
     plot(x=0:1,y=0:1,type="n",axes=FALSE,ann=FALSE)
     par(mai=c(0,0,0,0))
     if (is.null(shading))
      legend(0.02,1,levels(class.factor),mypalette[1:nclasses],cex=1.2)
     else
      legend(0.02,1,levels(class.factor),mypalette[1:nclasses],density=(1:nclasses) * shading,angle=rep(45,nclasses),cex=1.2)
     par(mai=mai.save)
    }
}

