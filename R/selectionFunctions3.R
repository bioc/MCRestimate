


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


varSel.svm.rfe <- function(sample.gene.matrix, classfactor,theParameter=NULL, ...)
{if(is.null(theParameter)){
    require(rfe)
    ## these functions are masking three buggy functions of the package rfe.
    assign("crossval",function (x, y, theta.fit, theta.predict, ..., nfold = n,nbmodels=1){
     call <- match.call()
     x <- as.matrix(x)
     n <- length(y)
     nfold <- trunc(nfold)
     if (nfold < 2) {
       stop("nfold should be greater than or equal to 2")
     }
     if (nfold > n) {
       stop("nfold should be less than or equal to the number of observations")
     }
     if (nfold == n) {
       groups <- 1:n
       leave.out <- 1
     }
     if (nfold < n) {
       leave.out <- trunc(n/nfold)
       o <- sample(1:n)
       groups <- vector("list", nfold)
       for (j in 1:(nfold - 1)) {
         jj <- (1 + (j - 1) * leave.out)
         groups[[j]] <- (o[jj:(jj + leave.out - 1)])
       }
       groups[[nfold]] <- o[(1 + (nfold - 1) * leave.out):n]
     }
     u <- vector("list", nfold)
     cv.fit <- data.frame(matrix(y, n,nbmodels))
     for (j in 1:nfold) {
       x=data.frame(x)
       u[[j]] <- theta.fit(x[-groups[[j]], ], y[-groups[[j]]], ...)
       cv.fit[groups[[j]],] <- theta.predict(u[[j]], x[groups[[j]], ])
     }
    #if (leave.out == 1)
    #groups <- NULL
    return( list(model.fit=u, cv.fit = cv.fit, nfold = nfold, leave.out = leave.out,groups = groups, call = call))
   },envir=.GlobalEnv)
           

   #
   #
    
   assign("rfe.predict",function(fit,x){
    nbmodels<-length(fit$TestedModels)
    d=as.matrix(x)
    if (dim(d)[2]==1){d=t(d)}
    yhat<-data.frame(matrix(0,dim(d)[1],nbmodels))
    for (j in 1:nbmodels){
    yhat[,j]<-predict((fit$TestedModels[[j]])$model,x[,(fit$TestedModels[[j]])$modelFeatures])}
    return(yhat)
   },envir=.GlobalEnv)

  assign("orderFeatures", function (x, y, FeaturesFromBest2Worst, nbFeatures, speed){
    p <- dim(x)[2]
    if (nbFeatures == p) {
        FeaturesToTest <- (1:p)
    }
    else {
        FeaturesToTest <- FeaturesFromBest2Worst[1:nbFeatures]
    }
    model <- svm(x[, FeaturesToTest], y,kernel="linear")
    sumabsw <- linear.feasel(model)
    FeatureOrder <- order(sumabsw, decreasing = TRUE)
    FeaturesFromBest2Worst[1:nbFeatures] <- FeaturesToTest[FeatureOrder]
    if (speed == "low") {
        nbFeatures <- nbFeatures - 1
    }
    else {
        if (nbFeatures == 2^(floor(log2(nbFeatures)))) {
            nbFeatures <- max(1, nbFeatures/2)
        }
        else {
            nbFeatures <- 2^(floor(log2(nbFeatures)))
        }
    }
    return(list(nbFeatures = nbFeatures, FeaturesFromBest2Worst =FeaturesFromBest2Worst,model = model, modelFeatures = FeaturesToTest))
    },envir=.GlobalEnv)
    
    x.t=t(sample.gene.matrix)
    cv.err = rfe.cv(x.t, classfactor)
    best.err = min(cv.err$error.cv)
    pos.best.err = match(best.err, cv.err$error.cv)
    no.features = cv.err$nbFeatures[pos.best.err]
    fit = rfe.fit(x.t, classfactor)
    index.sel.feat=fit$Flist[1:no.features]
    theParameter=rep(TRUE,nrow(sample.gene.matrix))
    theParameter[index.sel.feat]=FALSE
  }
    if (sum(!theParameter)==0) stop("The preprocessing method 'varSel.svm.rfe' chose no parameter.")
    train.matrix <- sample.gene.matrix[!theParameter,,drop=FALSE]
    return(list(matrix=train.matrix, parameter=theParameter))
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


