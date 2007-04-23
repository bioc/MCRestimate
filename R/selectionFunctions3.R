


####
#new methods for variable selection and NA replacement
#general interface:
#x: a matrix in which the rows corresponds to genes and  the colums corresponds to samples
#y: classfactor
#general function value
#a list with two named components:
#matrix: modified matrix x
#parameter: a logical vector with as many components as genes in the input matrix x, indicating by TRUE which genes have been removed for generation of the output matrix
  
#every gene that has a NA value in more than NAthreshold % of all samples is removed
varSel.removeManyNA <- function (sample.gene.matrix, classfactor, theParameter=NULL, NAthreshold=0.25,...)
{if(is.null(theParameter)){
   theParameter <- rep(FALSE, nrow(sample.gene.matrix))
   sum.na       <- sum(is.na(sample.gene.matrix))
   if (sum.na != 0) 
     theParameter <- select.NA.elements(sample.gene.matrix, NAthreshold = NAthreshold, byRow=TRUE)
   }
 x.removedNA <- sample.gene.matrix[!theParameter,,drop=FALSE]
 return(list(matrix=x.removedNA, parameter=theParameter))
}

varSel.impute.NA <- function(sample.gene.matrix ,classfactor, theParameter=NULL,...) 
{ if (is.null(theParameter))
   theParameter  <- apply(sample.gene.matrix, 1, median,na.rm=TRUE)
  cluster.gene.matrix <- replace.NA(sample.gene.matrix,theParameter,byRow=TRUE)  
  return(list(matrix=cluster.gene.matrix,parameter=theParameter))
}


replace.NA<- function(x, replacement , byRow = TRUE){
if (byRow==TRUE){
  norows=nrow(x)
  lengthy=length(replacement)
  if (norows==lengthy) {
    for (i in 1:norows) {
     x[i,is.na(x[i,])]=replacement[i]
    }
  } else {
print("Error: length of replacement vector does not match number of rows")
  }
  }
else {
nocol=ncol(x)
lengthy=length(replacement)
if (nocol==lengthy) {
for (i in 1:nocol) {
x[is.na(x[,i]),i]=replacement[i]
  }
  }
else {
print("Error: length of replacement vector does not match number of columns")
  }
  }
return(x)
}


############
#helper functions
############


##gibt eine vektor mit TRUE/FALSE zurück!! TRUE entspricht wohl, die zeile zu streichen

select.NA.elements <-
function(x, NAthreshold, byRow = TRUE){
  if (byRow){
    norow=nrow(x)
    nocol=ncol(x)
    y=NULL
    for (i in 1:norow){
      noNArow=sum(is.na(x[i,]))
      threshrow=noNArow/nocol
      if (threshrow>=NAthreshold){
        y=c(y,TRUE)
      }
      else {
        y=c(y,FALSE)
      }
    }
  }
  else {
    nocol=ncol(x)
    norow=nrow(x)
    y=NULL
    for (i in 1:nocol){
      noNAcol=sum(is.na(x[,i]))
      threshcol=noNAcol/norow
      if (threshcol>=NAthreshold){
        y=c(y,TRUE)
      }
    else {
      y=c(y,FALSE)
    }
   }
  }
return(y)
}


