identity <- function(sample.gene.matrix,classfactor,...) list(matrix=sample.gene.matrix,parameter=0)


#######################################
## 500 Gene with highest t-statistic ##
#######################################


# exprSet version

varSel.highest.t.stat <- function(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=500,...)
  {if (nlevels(classfactor)>2)stop("The gene reduction function 'varSel.highest.t.stat' only works with two classes. Please choose another gene reduction function or use a data set with only two classes.")
   if(is.null(theParameter))
 { require(arrayMagic)
   tscores                <- rowttests(sample.gene.matrix,classfactor)$statistic
   selection              <- order(abs(tscores),decreasing=TRUE)[1:var.numbers]
  theParameter            <- rep(TRUE,nrow(sample.gene.matrix))
  theParameter[selection] <- FALSE
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}


# eSRG version

varSel.highest.t.stat.eSRG <- function(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=500,...)
  {if (nlevels(classfactor)>2)stop("The gene reduction function 'varSel.highest.t.stat' only works with two classes. Please choose another gene reduction function or use a data set with only two classes.")
   if(is.null(theParameter))
 { require(arrayMagic)
   m                     <- nrow(sample.gene.matrix)/2
   new.matr              <- sample.gene.matrix[1:m,]-sample.gene.matrix[(m+1):(2*m),]
   tscores               <- rowttests(new.matr,classfactor)$statistic
   selection             <- order(abs(tscores),decreasing=TRUE)[1:var.numbers]
   bad.values            <- rep(TRUE,nrow(sample.gene.matrix))
   bad.values[selection] <- FALSE
   theParameter            <- c(bad.values,bad.values)
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}


#################################
## Genes with highest variance ##
#################################

# exprSet version

varSel.highest.var <- function(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=2000,...)
 {if(is.null(theParameter))
 {gene.sd   <- apply (sample.gene.matrix,1,var)
  selection <- order(gene.sd,decreasing=TRUE)[1:var.numbers]                   
  theParameter            <- rep(TRUE,nrow(sample.gene.matrix))
  theParameter[selection] <- FALSE
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}

# eSRG version

varSel.highest.var.eSRG <- function(sample.gene.matrix,classfactor,theParameter=NULL,var.numbers=2000,...)
 {if(is.null(theParameter))
 {m                     <- nrow(sample.gene.matrix)/2
  new.matr              <- sample.gene.matrix[1:m,]-sample.gene.matrix[(m+1):(2*m),]
  gene.sd               <- apply (new.matr,1,var)
  selection             <- order(gene.sd,decreasing=TRUE)[1:var.numbers]                   
  bad.values            <- rep(TRUE,m)
  bad.values[selection] <- FALSE
  theParameter <- c(bad.values,bad.values)
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}

###################################################################################

varSel.green.int.max.eSRG <- function(sample.gene.matrix,classfactor,theParameter=NULL, lambda=0.5,...)
 {if(is.null(theParameter))
 { m <- nrow(sample.gene.matrix)/2
   sset <- 1:m  
   gene.max   <- apply (sample.gene.matrix[sset,],1,max)
   a <- quantile(gene.max,lambda)
   bad.values <- gene.max < a
   theParameter <- c(bad.values,bad.values)
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}



varSel.green.int.sec.eSRG <- function(sample.gene.matrix,classfactor,theParameter=NULL, lambda=0.5,...)
 {if(is.null(theParameter))
 { m <- nrow(sample.gene.matrix)/2
   sset <- 1:m  
   gene.max   <- apply (sample.gene.matrix[sset,],1,function(x){sort(x,decreasing=TRUE)[2]})
   bad.values <- gene.max <  quantile(gene.max,lambda)
   theParameter <- c(bad.values,bad.values)
  }
  train.matrix   <- sample.gene.matrix[!theParameter,,drop=FALSE]
  return(list(matrix=train.matrix,parameter=theParameter))
}



varSel.AUC <- function (sample.gene.matrix, classfactor, theParameter = NULL,var.numbers=200,...)
{   require(ROC)
  if (nlevels(classfactor)>2)stop("The gene reduction function 'varSel.AUC' only works with two classes. Please choose another gene reduction function or use a data set with only two classes.")
    if (is.null(theParameter)) {
      newfactor <- as.factor(as.numeric(classfactor)-1)
      AUC.rfc<-function(gene.exprs){ return(AUC(rocdemo.sca(newfactor,gene.exprs,dxrule.sca,caseLabel=" ",markerLabel=" ")))}
   gene.AUC <- apply(sample.gene.matrix, 1, AUC.rfc)
   AUC.over.under.res<-as.vector(apply(cbind(gene.AUC,1-gene.AUC),1,max))
   theParameter <- !( AUC.over.under.res) %in% sort(AUC.over.under.res,decreasing=TRUE)[1:var.numbers]
    }
    train.matrix <- sample.gene.matrix[!theParameter, , drop = FALSE]
    return(list(matrix = train.matrix, parameter = theParameter))
}



cluster.kmeans.mean<- function(sample.gene.matrix ,classfactor, theParameter=NULL,number.clusters=500,...) 
{ if (is.null(theParameter)) theParameter <- kmeans(sample.gene.matrix,number.clusters,iter.max=50)
  number.of.cluster <- length(theParameter$size)
  cluster.gene.matrix <- matrix(NA,nrow=number.of.cluster, ncol=ncol(sample.gene.matrix))
  for (i in 1:number.of.cluster)
   cluster.gene.matrix[i,] <- colMeans(sample.gene.matrix[theParameter$cluster==i,,drop=FALSE])
  
  return(list(matrix=cluster.gene.matrix,parameter=theParameter))
}
