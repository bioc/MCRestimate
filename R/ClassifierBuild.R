ClassifierBuild <- function(eset,
                        class.column,
                        reference.class=NULL,
                        classification.fun,
                        variableSel.fun ="identity",
                        cluster.fun ="identity",
                        poss.parameters=list(),
                        cross.inner=10,
                        rand=123,
                        information=TRUE,
                        thePreprocessingMethods=c(variableSel.fun,cluster.fun)) UseMethod("ClassifierBuild")

ClassifierBuild.default <- function(eset,
                        class.column,
                        reference.class=NULL,
                        classification.fun,
                        variableSel.fun ="identity",
                        cluster.fun ="identity",
                        poss.parameters=list(),        
                        cross.inner=10,
                        rand=123,
                        information=TRUE,
                        thePreprocessingMethods=c(variableSel.fun,cluster.fun))

{ 
	# check the classification functions
  if(!(length(classification.fun)==1 & is.function(get(classification.fun))))
    stop("'classification.fun' must be an existing classification function.")

	# check the preprocessing methods
  for (k in 1:length(thePreprocessingMethods))
    if(!(is.function(get(thePreprocessingMethods[k]))))
     stop(paste("Preprocessing method",k,"must be an existing gene reduction function."))

	# check the parameters of classification function
  para.names.cla.fct <- setdiff(names(formals(classification.fun)),c("sample.gene.matrix","theParameter","classfactor","...","x","y"))
  para.names.sel.fct <- c()
	
  for (k in 1:length(thePreprocessingMethods))
  	para.names.sel.fct <- c(para.names.sel.fct,setdiff(names(formals(thePreprocessingMethods[k])),c("sample.gene.matrix","theParameter","classfactor","...","x","y")))
  
  if( !all (names(poss.parameters) %in% unique(c(para.names.sel.fct,para.names.cla.fct))))
    warning("There is at least one parameter in the parameter list that does not correspond to a parameter of the classification, cluster or selection function.")

  if ("theParameter" %in% names(poss.parameters))
    stop("You have to change the parameter because 'theParameter','cluster.check', and 'throw.away' are not allowed as parameter names.")

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
   { m <- t(x)
     for (i in 1:length(thePreprocessingMethods))
       { gene.reduce             <- get(thePreprocessingMethods[i])(m,label,theParameter=PreprocessingSlots[[i]],,...)
         m                       <- gene.reduce$matrix
         PreprocessingSlots[[i]] <- gene.reduce$parameter
       }
    return(list(transf.matrix=t(m),ThePreSlots=PreprocessingSlots))
   }
 the.function.for.classification <- function(x,y,...)
          {preprocessing.step <- preprocessing(x,y,...)
           tf.matrix          <- preprocessing.step$transf.matrix
           PreproParameter    <- preprocessing.step$ThePreSlots
           aaa                <- get(classification.fun)(tf.matrix,y,...)
           predict.function   <- aaa$predict
           assign("model.information",aaa$info,envir=we)

           rfct <- function(test.matrix)
           {testmatrix <- preprocessing(test.matrix,y,PreprocessingSlots=PreproParameter)$transf.matrix
            result     <- predict.function(testmatrix)
            return(result)
           }
          return(rfct)
          }
 
  predicted <- function(model,test) return(model(test))


 ###### the random generator is set
 if(! is.null(rand))
  set.seed(rand)


 the.vector.of.all.parameters <- vector(length=0, mode="list")
 eset.matrix                  <- exprs(eset)  # row=genes, col=samples
 Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model






 train.matrix <- t(eset.matrix) #row=samples,col=genes
 train.factor <- class.column.factor


# now the 'best' parameter are calculated if this is nessessary ( if there are choices)
     if (! (all(sapply(poss.parameters,length) %in% 1))) {
    result <- tune(the.function.for.classification,
	           train.matrix,
	           train.factor,
	           ranges=poss.parameters,
	           predict.func=predicted,
	           tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
    parameter.list <- result$best.parameter
    the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
 }else{
    parameter.list <- as.data.frame(poss.parameters)
 }
  
  predict.function <- do.call(the.function.for.classification, c(list(x=train.matrix,y=train.factor),parameter.list))
     

  classificationFunction <- function (eset){
    new.matr <- t(exprs(eset))
    if (ncol(train.matrix) != ncol(new.matr)) stop("The number of genes in the eset does not match the number of genes that are needed for the classifier.")
    return(predict.function(new.matr))
  }
 if(information)
  Model.information.List[[1]] <- with(we,model.information)



  # transforming the vector of parameters
  
  if (length(the.vector.of.all.parameters)!=0)
  {parameter.vector <- unlist (the.vector.of.all.parameters)
  the.vector.of.all.parameters <- tapply(as.vector(parameter.vector),factor(names(parameter.vector)),c)
  }


  
  rv <- list(classifier.for.matrix=predict.function,
             classifier.for.exprSet = classificationFunction,
             classes=class.column.factor,
             parameter=the.vector.of.all.parameters,
             thePreprocessingMethods=thePreprocessingMethods,
             class.method=classification.fun,
             cross.inner=cross.inner,
             information=Model.information.List)


#  environment(tune)  <-  original.tune.environment

  return(rv)
}

## Version for exprSetRG (arrayMagic)
#
# ClassifierBuild.exprSetRG<- function(eset,
#                         class.column,
#                         reference.class=NULL,
#                         classification.fun,
#                         variableSel.fun ="identity",
#                         cluster.fun ="identity",
#                         poss.parameters=list(), 
#                         cross.inner=10,
#                         rand=123,
#                         information=TRUE,
#                         thePreprocessingMethods=c(variableSel.fun,cluster.fun))
# 
# { require(arrayMagic)
#   if(!(length(classification.fun)==1 & is.function(get(classification.fun))))
#     stop("'classification.fun' must be an existing classification function.")
#   for (k in 1:length(thePreprocessingMethods))
#     if(!(is.function(get(thePreprocessingMethods[k]))))
#      stop(paste("Preprocessing method",k,"must be an existing gene reduction function."))
# 
#   para.names.cla.fct <- setdiff(names(formals(classification.fun)),c("sample.gene.matrix","theParameter","classfactor","...","x","y"))
#   para.names.sel.fct <- c()
#   for (k in 1:length(thePreprocessingMethods))
#     para.names.sel.fct <- c(para.names.sel.fct,setdiff(names(formals(thePreprocessingMethods[k])),c("sample.gene.matrix","theParameter","classfactor","...","x","y")))
# 
#   if( !all (names(poss.parameters) %in% unique(c(para.names.sel.fct,para.names.cla.fct))))
#     warning("There is at least one parameter in the parameter list that does not correspond to a parameter of the classification, cluster or selection function.")
# 
#   if ("theParameter" %in% names(poss.parameters))
#     stop("You have to change the parameter because 'theParameter','cluster.check' and 'throw.away' are not allowed as parameter names.")
# 
# 
# 
#   class.column.factor <- class.factor.format(getExprSetGreenMinusRed(eset),class.column,reference.class)
#   levels.class        <- levels(class.column.factor)
#   nlevels.class       <- nlevels(class.column.factor)
#   nn                  <- length(class.column.factor)
# 
# 
#   we                        <- new.env()
#   #original.tune.environment <- environment(tune)
#   #environment(tune)         <- we
# 
#   ## x is a matrix with rows=samples, columns=genes (for tune)
#   ## m is a matrix with columns=samples, rows=genes (as usual with exprSets)
#   tIp <- vector(mode="list",length=length(thePreprocessingMethods))
# 
# 
# 
# 
#   preprocessing <- function(x,label,PreprocessingSlots=tIp,...)
#    { m <- t(x)
#      for (i in 1:length(thePreprocessingMethods))
#        { gene.reduce             <- get(thePreprocessingMethods[i])(m,label,theParameter=PreprocessingSlots[[i]],,...)
#          m                       <- gene.reduce$matrix
#          PreprocessingSlots[[i]] <- gene.reduce$parameter
#        }
#     return(list(transf.matrix=t(m),ThePreSlots=PreprocessingSlots))
#    }
# 
# 
#   the.function.for.classification <-   function(x,y,...)
#           {preprocessing.step <- preprocessing(x,y,...)
#            tf.matrix          <- preprocessing.step$transf.matrix
#            a                  <- ncol(tf.matrix)
#            tf.matrix          <- tf.matrix[,1:(a/2)] - tf.matrix[,(a/2+1):a]  # now we have the log ratios for further things
#            PreproParameter    <- preprocessing.step$ThePreSlots
#            aaa                <- get(classification.fun)(tf.matrix,y,...)
#            predict.function   <- aaa$predict
#            assign("model.information",aaa$info,envir=we)
# 
#            rfct <- function(test.matrix)
#            {testmatrix <- preprocessing(test.matrix,y,PreprocessingSlots=PreproParameter)$transf.matrix
#             a          <- ncol(testmatrix)
#             testmatrix <- testmatrix[,1:(a/2)] - testmatrix[,(a/2+1):a]
#             result     <- predict.function(testmatrix)
#             return(result)
#            }
#           return(rfct)
#           }
# 
# 
# 
# 
# 
# 
#   predicted <- function(model,test) return(model(test))
# 
# 
#  ###### the random generator is set
#  if(! is.null(rand))
#   set.seed(rand)
# 
# 
# 
#  the.vector.of.all.parameters <- vector(length=0, mode="list")
#  eset.matrix              <- rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset)))  # row=genes, col=samples genes are two times
#  Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model
# 
# 
#  train.matrix <- t(eset.matrix) #row=samples,col=genes
#  train.factor <- class.column.factor
# 
#   if (! (all(sapply(poss.parameters,length) %in% 1))){
#      result<-tune(the.function.for.classification,
#             train.matrix,
#             train.factor,
#             ranges=poss.parameters,
#             predict.func=predicted,
#             tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
#      parameter.list <- result$best.parameter
#      the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
#  } else {
#      parameter.list <- as.data.frame(poss.parameters)
#      }
#  predict.function <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))
# 
# 
#  classificationFunction <- function (eset){
#    new.matr <- t(rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset))))
#    if (ncol(train.matrix) != ncol(new.matr)) stop("The number of genes in the eset does not match the number of genes that are needed for the classifier.")
#    return(predict.function(new.matr))
#  }
# 
# 
#  if (length(the.vector.of.all.parameters)!=0)
#   {parameter.vector             <- unlist (the.vector.of.all.parameters)
#    the.vector.of.all.parameters <- tapply(as.vector(parameter.vector),factor(names(parameter.vector)),c)
#   }
# 
# 
# 
# ## all data are collected
# if(information)  
#  Model.information.List[[1]] <- with(we,model.information)
# 
# 
#   rv <- list(classifier=predict.function,
#              classifier.for.exprSetRG= classificationFunction,
#              classes=class.column.factor,
#              parameter=the.vector.of.all.parameters,
#              thePreprocessingMethods=thePreprocessingMethods,
#              class.method=classification.fun,
#              ross.inner=cross.inner,
#              information=Model.information.List)
# 
#   #environment(tune)  <-  original.tune.environment
# 
#   return(rv)
# }
