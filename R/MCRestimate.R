
MCRestimate <- function(eset,
                        class.column,
                        reference.class=NULL,
                        classification.fun,
                        variableSel.fun ="identity",
                        cluster.fun ="identity",
                        poss.parameters=list(),
                        cross.outer=10,
                        cross.repeat=3,
                        cross.inner=cross.outer,
                        plot.label=NULL,
                        rand=123,
                        stratify=FALSE,
                        information=TRUE,
                        block.column=NULL,
                        thePreprocessingMethods=c(variableSel.fun,cluster.fun)) UseMethod("MCRestimate")

MCRestimate.default <- function(eset,
                        class.column,
                        reference.class=NULL,
                        classification.fun,
                        variableSel.fun ="identity",
                        cluster.fun ="identity",
                        poss.parameters=list(), 
                        cross.outer=10,
                        cross.repeat=3,
                        cross.inner=cross.outer,
                        plot.label=NULL,
                        rand=123,
                        stratify=FALSE,
                        information=TRUE,
                        block.column=NULL,
                        thePreprocessingMethods=c(variableSel.fun,cluster.fun))

{ 
	if(!(length(classification.fun)==1 & is.function(get(classification.fun))))
    stop("'classification.fun' must be an existing classification function.")
  for (k in 1:length(thePreprocessingMethods))
    if(!(is.function(get(thePreprocessingMethods[k]))))
     stop(paste("Preprocessing method",k,"must be an existing gene reduction function."))

  para.names.cla.fct <- setdiff(names(formals(classification.fun)),c("sample.gene.matrix","theParameter","classfactor","...","x","y"))
  para.names.sel.fct <- c()
  for (k in 1:length(thePreprocessingMethods))
  para.names.sel.fct <- c(para.names.sel.fct,setdiff(names(formals(thePreprocessingMethods[k])),c("sample.gene.matrix","theParameter","classfactor","...","x","y")))

  if( !all (names(poss.parameters) %in% unique(c(para.names.sel.fct,para.names.cla.fct))))
    warning("There may be at least one parameter in the parameter list that does not correspond to a parameter of the classification or any preprocessing function.")

  if ("theParameter" %in% names(poss.parameters))
    stop("You have to change the parameter because 'theParameter','cluster.check', and 'throw.away' are not allowed as parameter names.")

	# the method class.factor.format returns the desired column from the pData of eset.
	# parameter ref.class is set if the user wishes to	compare one class against all other classes
  class.column.factor <- class.factor.format(eset,class.column,reference.class)
  levels.class        <- levels(class.column.factor)
  nlevels.class       <- nlevels(class.column.factor)
  nn                  <- length(class.column.factor)
  
  we                        <- new.env()
  #original.tune.environment <- environment(tune)
  #environment(tune)         <- we

  ## x is a matrix with rows=samples, columns=genes (for tune)
  ## m is a matrix with columns=samples, rows=genes (as usual with exprSets)
  tIp <- vector(mode="list",length=length(thePreprocessingMethods))

  preprocessing <- function(x,label,PreprocessingSlots=tIp,...)
  {
		m <- t(x)
		for (i in 1:length(thePreprocessingMethods))
		{
				 gene.reduce             <- get(thePreprocessingMethods[i])(m,label,theParameter=PreprocessingSlots[[i]],,...)
	       m                       <- gene.reduce$matrix
	       PreprocessingSlots[[i]] <- gene.reduce$parameter
		}
		return(list(transf.matrix=t(m),ThePreSlots=PreprocessingSlots))
	}
 

  the.function.for.classification<-  function(x,y,...)
	{
		preprocessing.step <- preprocessing(x,y,...)
		tf.matrix          <- preprocessing.step$transf.matrix
		PreproParameter    <- preprocessing.step$ThePreSlots
		aaa                <- get(classification.fun)(tf.matrix,y,...)
		predict.function   <- aaa$predict
		assign("model.information",aaa$info,envir=we)
			
		rfct <- function(test.matrix)
		{
			testmatrix <- preprocessing(test.matrix,y,PreprocessingSlots=PreproParameter)$transf.matrix
			result     <- predict.function(testmatrix)
			return(result)
		}
		return(rfct)
	}


  predicted <- function(model,test) return(model(test))

  my.balanced.folds <- function(class.column.factor, cross.outer){
    sampleOfFolds  <- get("balanced.folds",en=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
    permutated.cut <- rep(0,length(class.column.factor))
    for (sample in 1:cross.outer){
      permutated.cut[sampleOfFolds[[sample]]] <- sample
    }
    return(permutated.cut)
  }


  
 ###### the random generator is set
 if(! is.null(rand))
  set.seed(rand)

 #######creating the block of samples which should be in the test set


 if (cross.outer > nn | cross.outer < 2 ) stop(paste("You have to specify at least two different groups and not more than",nn))
 if (cross.repeat < 1)        stop("There must be at least one repetition of the cross validation procedure")

 the.cut <- cut(1:nn,cross.outer,1:cross.outer)

 ROWNAMES                     <- class.column.factor
 COLNAMES                     <- levels.class
 the.vector.of.all.parameters <- vector(length=0, mode="list")

 eset.matrix                  <- exprs(eset)  # row=genes, col=samples

  #assign("model.information",list(),envir=we)
  Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model

  if( ! is.null(block.column) ){
    if (!block.column %in% names(pData(eset))){
      stop("The value of 'block.column' should be the name of a column of pData(eset).")
    }
    block.factor <- as.factor(pData(eset)[,block.column])
    the.blocks <-  as.integer(block.factor)
    if( ! all(tapply(as.vector(class.column.factor),the.blocks, function(x) length(unique(x))))==1){
      stop(" inconsistent class assignment for the different blocks in block.column ")
    }
    class.column.factor.blocks <- class.column.factor[!duplicated(the.blocks)]
    nb <- length(unique(the.blocks))
    if( cross.outer > nb ){
      warning(paste("cross.outer: You have to specify not more than",nb," i.e. the number of different blocks in block.column; cross.outer is set to this number"))
      cross.outer <- nb
    }
    if( cross.outer == nb ){
      if( cross.repeat != 1 ){
        cross.repeat <- 1
        warning("cross.outer equals the number of blocks in block.column and hence cross.repeat is set to 1")     }
    }
    the.block.cut <-  cut(1:nb,cross.outer,1:cross.outer)
  }else{
    block.factor <- NULL
  }

  ## votal.matrix                <- matrix(0, ncol=nlevels.class, nrow=nn)
  AllVotes                     <- array (NA, dim=c(nn,nlevels.class, cross.repeat))
  permutated.cut.matrix        <- array (NA, dim=c(nn, cross.repeat))



  for (l in 1:cross.repeat){
    the.votes.per.cv             <- matrix(0, ncol=nlevels.class,nrow=nn)
    rownames(the.votes.per.cv)   <- class.column.factor
    colnames(the.votes.per.cv)   <- levels.class
    if( ! is.null(block.column) ){

      if( stratify ){
        permutated.block.cut <- my.balanced.folds(class.column.factor.blocks, cross.outer)
      }else{
        permutated.block.cut <- sample(the.block.cut,nb)
      }
      permutated.cut <- rep(NA, nn)
      for (sample in 1:cross.outer){
        blockSel <- which(permutated.block.cut == sample)
        blocks <- the.blocks[!duplicated(the.blocks)][blockSel]
        permutated.cut[the.blocks %in% blocks] <- sample
      }      
    }else{
      if (stratify){
        permutated.cut <- my.balanced.folds(class.column.factor, cross.outer)
      }else{
        permutated.cut <- sample(the.cut,nn)
      }
    }
        
  ## the start of the outer cross-validation ##
  #############################################

    for (sample in 1:cross.outer)
     {

   # dividing the set into a testset and a trainingsset.
     block <- permutated.cut==sample
         
     train.matrix <- t(eset.matrix[,!block,drop=FALSE]) #row=samples,col=genes
     test.matrix  <- t(eset.matrix[,block,drop=FALSE])
     train.factor <- class.column.factor[!block]


   # now the 'best' parameter are calculated if this is nessessary ( if there are choices)
     if (! (all(sapply(poss.parameters,length) %in% 1)))
      {#if(length(train.factor)%/%cross.inner <= 1) stop(paste("Please choose a value for 'cross.inner' with the following attribute:",length(train.factor),"/ 'cross.inner' > 2"))
        result <- tune(the.function.for.classification,
                         train.matrix,
		         train.factor,
                         ranges=poss.parameters,
                         predict.func=predicted,
                         tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
        parameter.list               <- result$best.parameter
        the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
        #predict.function <- result$best.model
       }
     else{
       parameter.list   <- as.data.frame(poss.parameters)
     }
       predict.function <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))


     
   # the votes for a class are calculated

    # predict.function         <- with(we,do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list)))
    if(information)
     Model.information.List[[(l-1)*cross.outer + sample]] <- with(we,model.information)
    
     pred.vector              <- predict.function(test.matrix)
     vote.matrix              <- t(sapply(1:length(pred.vector), function(j) as.numeric(levels.class==pred.vector[j])))
     colnames(vote.matrix)    <- levels.class
     the.votes.per.cv[block,] <- vote.matrix
   }
  permutated.cut.matrix[,l] <- permutated.cut
  AllVotes[,,l] <- the.votes.per.cv
  #votal.matrix <- votal.matrix + the.votes.per.cv
  }
  votal.matrix <- apply(AllVotes,c(1,2),sum)
  rownames(votal.matrix) <- ROWNAMES
  colnames(votal.matrix) <- COLNAMES
  votal.matrix           <- votal.matrix /cross.repeat
  
  # transforming the vector of parameters
  
  if (length(the.vector.of.all.parameters)!=0)
  {parameter.vector <-unlist(lapply(the.vector.of.all.parameters,as.character))
  the.vector.of.all.parameters <- tapply(as.vector(parameter.vector),factor(names(parameter.vector)),c)
  }
  # creating  the confusion table

  res <- whatiscorrect(votal.matrix)
  vote.table           <- table(rownames(votal.matrix), res$best.vote)
  new.table            <- matrix(0,
                                 ncol=nrow(vote.table),
                                 nrow=nrow(vote.table),
                                 dimnames=list(rownames(vote.table),rownames(vote.table)))
  
  new.table[,colnames(vote.table)] <- vote.table
                                      
  normed.table         <- new.table/rowSums(new.table)
  confusion            <- cbind(new.table, 1-diag(normed.table))
  colnames(confusion)  <- c(levels(class.column.factor), "class error")
 


## if plot.label is specified the sample.names are set
  
  if (!(is.null(plot.label)))
    if (length(plot.label)==1 ) sample.names <- pData(eset)[,plot.label]
    else                        sample.names <- plot.label
  else
    sample.names <- 1:nn

## all data are collected
  
  rv <- list(votes=votal.matrix,
                                classes=class.column.factor,
                                table=confusion,
                                correct.prediction=res$correct.prediction,
                                correct.class.vote=res$correct.class.vote,
                                parameter=the.vector.of.all.parameters,
                                class.method=classification.fun,
                                thePreprocessingMethods=thePreprocessingMethods,
                                cross.outer=cross.outer,
                                cross.inner=cross.inner,
                                cross.repeat=cross.repeat,
                                sample.names=sample.names,
                                information=Model.information.List,
                                stratify=stratify,
                                block.column=block.column,
                                block.factor=block.factor,
                                permutated.cut.matrix=permutated.cut.matrix,
                                indVotes=AllVotes)
                           
   
  class(rv) <- "MCRestimate"
  #environment(tune)  <-  original.tune.environment
  
  return(rv)
}


## Version for exprSetRG (arrayMagic)

#MCRestimate.exprSetRG<- function(eset,
#                        class.column,
#                        reference.class=NULL,
#                        classification.fun,
#                        variableSel.fun ="identity",
#                        cluster.fun ="identity",
#                        poss.parameters=list(),    
#                        cross.outer=10,
#                        cross.repeat=3,
#                        cross.inner=cross.outer,
#                        plot.label=NULL,
#                        rand=123,
#                        stratify=FALSE,
#                        information=TRUE,
#                        block.column=NULL,
#                        thePreprocessingMethods=c(variableSel.fun,cluster.fun))
#
#{ require(arrayMagic)
#  if(!(length(classification.fun)==1 & is.function(get(classification.fun))))
#    stop("'classification.fun' must be an existing classification function.")
#   for (k in 1:length(thePreprocessingMethods))
#    if(!(is.function(get(thePreprocessingMethods[k]))))
#     stop(paste("Preprocessing method",k,"must be an existing gene reduction function."))
#
#  para.names.cla.fct <- setdiff(names(formals(classification.fun)),c("sample.gene.matrix","theParameter","classfactor","...","x","y"))
#  para.names.sel.fct <- c()
#  for (k in 1:length(thePreprocessingMethods))
#    para.names.sel.fct <- c(para.names.sel.fct,setdiff(names(formals(thePreprocessingMethods[k])),c("sample.gene.matrix","theParameter","classfactor","...","x","y")))
#
#
#  if( !all (names(poss.parameters) %in% unique(c(para.names.sel.fct,para.names.cla.fct))))
#    warning("There may be at least one parameter in the parameter list that does not correspond to a parameter of the classification, cluster or selection function.")
#
#  if ("theParameter" %in% names(poss.parameters))
#    stop("You have to change the parameter because 'theParameter' is not allowed as parameter names.")
#  
#  
#  GreenMinusRedSet    <- getExprSetGreenMinusRed(eset) # it is not important if this is green minus red or red minus green because we are only interested in the phenoData.
#  class.column.factor <- class.factor.format(GreenMinusRedSet,class.column,reference.class)
#  levels.class        <- levels(class.column.factor)
#  nlevels.class       <- nlevels(class.column.factor)
#  nn                  <- length(class.column.factor)
#
#  
#  we                        <- new.env()
#  #original.tune.environment <- environment(tune)
#  #environment(tune)         <- we
#
#
#  
#  ## x is a matrix with rows=samples, columns=genes (for tune)
#  ## m is a matrix with columns=samples, rows=genes (as usual with exprSets)
#  tIp <- vector(mode="list",length=length(thePreprocessingMethods))
#
#  preprocessing <- function(x,label,PreprocessingSlots=tIp,...)
#   { m <- t(x)
#     for (i in 1:length(thePreprocessingMethods))
#       { gene.reduce             <- get(thePreprocessingMethods[i])(m,label,theParameter=PreprocessingSlots[[i]],,...)
#         m                       <- gene.reduce$matrix
#         PreprocessingSlots[[i]] <- gene.reduce$parameter
#       }
#    return(list(transf.matrix=t(m),ThePreSlots=PreprocessingSlots))
#   }
#
#
#  the.function.for.classification<-  function(x,y,...)
#          {preprocessing.step <- preprocessing(x,y,...)
#           tf.matrix          <- preprocessing.step$transf.matrix
#           PreproParameter    <- preprocessing.step$ThePreSlots
#           aaa                <- get(classification.fun)(tf.matrix,y,...)
#           predict.function   <- aaa$predict
#           assign("model.information",aaa$info,envir=we)
#
#           rfct <- function(test.matrix)
#           {testmatrix <- preprocessing(test.matrix,y,PreprocessingSlots=PreproParameter)$transf.matrix
#            result     <- predict.function(testmatrix)
#            return(result)
#           }
#          return(rfct)
#          }
#
#
#  predicted <- function(model,test) return(model(test))
#
#  my.balanced.folds <- function(class.column.factor, cross.outer){
#    sampleOfFolds  <- get("balanced.folds",en=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
#    permutated.cut <- rep(0,length(class.column.factor))
#    for (sample in 1:cross.outer){
#      permutated.cut[sampleOfFolds[[sample]]] <- sample
#    }
#    return(permutated.cut)
#  }
#
#
# ###### the random generator is set
# if(! is.null(rand))
#  set.seed(rand)
# 
# #######creating the block of samples which should be in the test set
# 
#
# if (cross.outer > nn | cross.outer < 2 ) stop(paste("You have to specify at least two different groups and not more than",nn))
# if (cross.repeat < 1)        stop("There must be at least one repetition of the cross validation procedure")
#
# the.cut <- cut(1:nn,cross.outer,1:cross.outer)
#
# ROWNAMES                     <- class.column.factor
# COLNAMES                     <- levels.class
# the.vector.of.all.parameters <- vector(length=0, mode="list")
#
# eset.matrix              <- rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset)))  # row=genes, col=samples genes are two times
#
# #assign("model.information",list(),envir=we) 
# Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model
#  
#  if( ! is.null(block.column) ){
#    if (!block.column %in% names(pData(eset))){
#      stop("The value of 'block.column' should be the name of a column of pData(eset).")
#    }
#    block.factor <- as.factor(pData(eset)[,block.column])
#    the.blocks <-  as.integer(block.factor)
#    if( ! all(tapply(as.vector(class.column.factor),the.blocks, function(x) length(unique(x))))==1){
#      stop(" inconsistent class assignment for the different blocks in block.column ")
#    }
#    class.column.factor.blocks <- class.column.factor[!duplicated(the.blocks)]
#    nb <- length(unique(the.blocks))
#    if( cross.outer > nb ){
#      warning(paste("cross.outer: You have to specify not more than",nb," i.e. the number of different blocks in block.column; cross.outer is set to this number"))
#      cross.outer <- nb
#    }
#    if( cross.outer == nb ){
#      if( cross.repeat != 1 ){
#        cross.repeat <- 1
#        warning("cross.outer equals the number of blocks in block.column and hence cross.repeat is set to 1")     }
#    }
#    the.block.cut <-  cut(1:nb,cross.outer,1:cross.outer)
#  }else{
#    block.factor <- NULL
#  }
#
#  ## votal.matrix                <- matrix(0, ncol=nlevels.class, nrow=nn)
#  AllVotes                     <- array (NA, dim=c(nn,nlevels.class, cross.repeat))
#  permutated.cut.matrix        <- array (NA, dim=c(nn, cross.repeat))
#
#
#  
#  for (l in 1:cross.repeat)
#  { the.votes.per.cv             <- matrix(NA, ncol=nlevels.class,nrow=nn)
#    rownames(the.votes.per.cv)   <- class.column.factor
#    colnames(the.votes.per.cv)   <- levels.class
#    if( ! is.null(block.column) ){
#
#      if( stratify ){
#        permutated.block.cut <- my.balanced.folds(class.column.factor.blocks, cross.outer)
#      }else{
#        permutated.block.cut <- sample(the.block.cut,nb)
#      }
#      permutated.cut <- rep(NA, nn)
#      for (sample in 1:cross.outer){
#        blockSel <- which(permutated.block.cut == sample)
#        blocks <- the.blocks[!duplicated(the.blocks)][blockSel]
#        permutated.cut[the.blocks %in% blocks] <- sample
#      }      
#    }else{
#      if (stratify){
#        permutated.cut <- my.balanced.folds(class.column.factor, cross.outer)
#      }else{
#        permutated.cut <- sample(the.cut,nn)
#      }
#    }
#      
#  ## the start of the outer cross-validation ##
#  #############################################
#    
#    for (sample in 1:cross.outer)
#     {
#
#   # dividing the set into a testset and a trainingsset.
#     block <- permutated.cut==sample
#
#     train.matrix <- t(eset.matrix[,!block,drop=FALSE]) #row=samples,col=genes
#     test.matrix  <- t(eset.matrix[,block,drop=FALSE])
#     train.factor <- class.column.factor[!block]
#
#
#   # now the 'best' parameter are calculated if this is nessessary ( if there are choices)
#     
#     if (! (all(sapply(poss.parameters,length) %in% 1)))
#      {#if(length(train.factor)%/%cross.inner <= 1) stop(paste("Please choose a value for 'cross.inner' with the following attribute:",length(train.factor),"/ 'cross.inner' > 2"))
#        result <- tune(the.function.for.classification,
#                         train.matrix,
#		         train.factor,
#                         ranges=poss.parameters,
#                         predict.func=predicted,
#                         tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
#       parameter.list <- result$best.parameter
#       the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
#       #predict.function <- result$best.model
#     }else{
#       parameter.list <- as.data.frame(poss.parameters)
#     }
#    predict.function         <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))
#
#
#   # the votes for a class are calculated
#   if(information)
#     Model.information.List[[(l-1)*cross.outer + sample]] <- with(we,model.information)
#     
#     pred.vector              <- predict.function(test.matrix)
#     vote.matrix              <- t(sapply(1:length(pred.vector), function(j) as.numeric(levels.class==pred.vector[j])))
#     colnames(vote.matrix)    <- levels.class
#     the.votes.per.cv[block,] <- vote.matrix
#
#   }
#  permutated.cut.matrix[,l] <- permutated.cut
#  AllVotes[,,l] <- the.votes.per.cv
#  #votal.matrix <- votal.matrix + the.votes.per.cv
#  }
#  votal.matrix <- apply(AllVotes,c(1,2),sum)
#  rownames(votal.matrix) <- ROWNAMES
#  colnames(votal.matrix) <- COLNAMES
#  votal.matrix <- votal.matrix /cross.repeat
#
#  # transforming the vector of parameters
#
#  if (length(the.vector.of.all.parameters)!=0)
#  {parameter.vector <-unlist(lapply(the.vector.of.all.parameters,as.character))
#  the.vector.of.all.parameters <- tapply(as.vector(parameter.vector),factor(names(parameter.vector)),c)
#  }
#  
#  # creating  the confusion table
#
#  res <- whatiscorrect(votal.matrix)
#  
#  
#  vote.table           <- table(rownames(votal.matrix), res$best.vote)
#  new.table            <- matrix(0,
#                                 ncol=nrow(vote.table),
#                                 nrow=nrow(vote.table),
#                                 dimnames=list(rownames(vote.table),rownames(vote.table)))
#  
#  new.table[,colnames(vote.table)] <- vote.table
#                                      
#  normed.table         <- new.table/rowSums(new.table)
#  confusion            <- cbind(new.table, 1-diag(normed.table))
#  colnames(confusion)  <- c(levels(class.column.factor), "class error")
# 
#
#
### if plot.label is specified 
#  
#  if (!(is.null(plot.label)))
#    if (length(plot.label)==1 ) sample.names <- pData(GreenMinusRedSet)[,plot.label]
#    else                        sample.names <- plot.label
#  else
#    sample.names <- 1:nn
#
### all data are collected
#  
#  rv <- list(votes=votal.matrix,
#                                classes=class.column.factor,
#                                table=confusion,
#                                correct.prediction=res$correct.prediction,
#                                correct.class.vote=res$correct.class.vote,
#                                parameter=the.vector.of.all.parameters,
#                                class.method=classification.fun,
#                                thePreprocessingMethods=thePreprocessingMethods,
#                                cross.outer=cross.outer,
#                                cross.inner=cross.inner,
#                                cross.repeat=cross.repeat,
#                                sample.names=sample.names,
#                                information=Model.information.List,
#                                stratify=stratify,
#                                block.column=block.column,
#                                block.factor=block.factor,
#                                permutated.cut.matrix=permutated.cut.matrix,
#                                indVotes=AllVotes)
#                           
#   
#  class(rv) <- "MCRestimate"
# # environment(tune)  <-  original.tune.environment
#  
#  return(rv)
#}
