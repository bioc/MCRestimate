MCRestimateMerge <- function(MCRestimateList){
## check if things are identical that should be identical
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

## Merging

## Parameter
thePara <- lapply(MCRestimateList, function(x) x$parameter)
stopifnot(all(sapply(thePara, function(x)identical(names(x),names(thePara[[1]])))))

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
colnames(confusion)  <- c(levels(MCRestimateList[[1]]$classes), "class error")

result <- MCRestimateList[[1]]
result$votes              <- votal.matrix
result$cross.repeat       <- rescr
result$table              <- confusion
result$correct.prediction <- res$correct.prediction
result$correct.class.vote <- res$correct.class.vote
result$information        <- resIn
result$parameter          <- resPara
return(result)
}
