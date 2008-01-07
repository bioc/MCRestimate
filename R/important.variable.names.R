important.variable.names <- function(mcr,file="important_variables",listName=NULL,writeFile=TRUE,...){
  if( !is.null(listName) ){
    #gene.list <- lapply(mcr$information,function(x) x[[list(...)$listName]])
    gene.list <- lapply(mcr$information,function(x) x[[listName]])
  }else{
    gene.list <- lapply(mcr$information,function(x) x[,1])
  }
  xx <- as.matrix(table (unlist(gene.list)))
  gene.names <- rownames(xx)
  the.length <- sapply(gene.list,length)
  ff  <- table(the.length)
  
  ## making a plot with group size frequency
  if(writeFile ){
    pdf(file=paste(file,"_freq.pdf",sep=""))
    if (length(ff) < 11) {gg <- as.vector(ff)
                          names(gg) <- names(ff)
                          barplot(gg,...)}
    else hist(the.length,br=ceiling(length(the.length)/5),col="yellow",...)
    dev.off()
  }
    
  positions <- list()
  kleinste.gruppe <- c()
  frq <- c()
  for (i  in 1:length(gene.names)){
   positions[[i]]     <- unlist(lapply(gene.list, function (x) match(gene.names[i],x)))
   kleinste.gruppe[i]<- min(the.length[sapply(gene.list, function(x) gene.names[i] %in% x)])}
  aa <-   sapply(positions,function(x) round(median(x,na=TRUE),2))
  bb <-   sapply(positions,function(x) round(mean(x,na=TRUE),2))
  cc <-   sapply(positions,function(x) round(var(x, na=TRUE),2))
  zz <- cbind(xx,aa,bb,cc,kleinste.gruppe)
  r <-  order(-zz[,1],zz[,2])
  zz <- zz[r,]
  colnames(zz) <- c("No. of times","Median position" , "Mean position","Variance","smallest group")
  if( writeFile ){
    write.table(zz,file=paste(file,".txt",sep=""),sep="\t",col.names = NA)
    write.table(as.matrix(ff),file=paste(file,"_freq.txt",sep=""),sep="\t",col.names = NA)
  }
  return(zz)
}

