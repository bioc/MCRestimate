MCRestimateMerge <- function(MCRestimateList){

if( length(MCRestimateList) == 1 ){
  return(MCRestimateList[[1]])
}
  
#### check if things are identical that should be identical

aa <- lapply(MCRestimateList, function(x) x$classes)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$class.method)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$select.method)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$cluster.method)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$cross.outer)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$cross.inner)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) x$sample.names)
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) rownames(x$votes))
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
aa <- lapply(MCRestimateList, function(x) colnames(x$votes))
stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))

##
varNames <- c("stratify","block.column", "block.factor")
for(varName in varNames){
  if( varName %in% as.vector(sapply(MCRestimateList, names))){
    aa <- lapply(MCRestimateList, function(x) x[[varName]])
    stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
  }
}

#### Merging


## Parameter
thePara <- lapply(MCRestimateList, function(x) x$parameter)
stopifnot(all(sapply(thePara, function(x)identical(names(x),names(thePara[[1]])))))



## indVotes

theindVotes <- lapply(MCRestimateList, function(x) x$indVotes)

xDim        <- unique(sapply(theindVotes, function(x) dim(x)[1]))
stopifnot(length(xDim)==1)

yDim        <- unique(sapply(theindVotes, function(x) dim(x)[2]))
stopifnot(length(yDim)==1)

zDim        <- sapply(theindVotes, function(x) dim(x)[3])
newA <- array(NA, dim=c(xDim,yDim,sum(zDim)))

newA[,,1:zDim[1]] <- theindVotes[[1]]
for(j in 2:length(zDim)){
 Skip <-  sum(zDim[1:(j-1)])
 newA[,,(Skip + 1):(Skip + zDim[j])] <- theindVotes[[j]]
}

## Information
theInfo <- lapply(MCRestimateList, function(x) x$information)

## Votes
aa <- lapply(MCRestimateList, function(x) x$votes)
stopifnot(all(sapply(aa, function(x)identical(rownames(x),rownames(aa[[1]])))))

bb <- lapply(MCRestimateList, function(x) x$cross.repeat)
dd <- mapply(function(x,y) x*y,aa,bb,SIMPLIFY=FALSE)

resv    <- dd[[1]]
rescr   <- bb[[1]]
resIn   <- theInfo[[1]]
resPara <- thePara[[1]]

for( i in 2:length(MCRestimateList)){
 
 resv <- resv + dd[[i]]
 rescr<- rescr + bb[[i]]
 resIn   <- c(resIn  , theInfo[[i]])
 resPara <- mapply(function(x,y) c(x,y),resPara,thePara[[i]],SIMPLIFY=FALSE)
}
names(resPara) <- names(thePara[[1]])
votal.matrix <- resv /rescr
res          <- whatiscorrect(votal.matrix)
vote.table   <- table(rownames(votal.matrix), res$best.vote)
new.table    <- matrix(0,ncol=nrow(vote.table),
                       nrow=nrow(vote.table),
                       dimnames=list(rownames(vote.table),rownames(vote.table)))
new.table[,colnames(vote.table)] <- vote.table
normed.table         <- new.table/rowSums(new.table)
confusion            <- cbind(new.table, 1-diag(normed.table))
colnames(confusion)  <- c(colnames(new.table), "class error")

result <- MCRestimateList[[1]]
result$votes              <- votal.matrix
result$cross.repeat       <- rescr
result$table              <- confusion
result$correct.prediction <- res$correct.prediction
result$correct.class.vote <- res$correct.class.vote
result$information        <- resIn
result$parameter          <- resPara
result$indVotes           <- newA

varName <- "permutated.cut.matrix"
if( varName %in% as.vector(sapply(MCRestimateList, names))){
  aa <- lapply(MCRestimateList, function(x) dim(x[[varName]]) )
  stopifnot(all(sapply(aa, function(x) identical(x,aa[[1]]))))
  new.permutated.cut.matrix <- do.call("cbind", args=lapply(MCRestimateList, function(x) x[[varName]]))
  result[[varName]] <- new.permutated.cut.matrix
}

return(result)
}
