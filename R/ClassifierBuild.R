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

{ if(!(length(classification.fun)==1 & is.function(get(classification.fun))))
    stop("'classification.fun' must be an existing classification function.")
  for (k in 1:length(thePreprocessingMethods))
    if(!(is.function(get(thePreprocessingMethods[k]))))
     stop(paste("Preprocessing method",k,"must be an existing gene reduction function."))


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



Mytune<-function (method, train.x, train.y = NULL, data = list(), validation.x = NULL,
    validation.y = NULL, ranges = NULL, predict.func = predict,
    tunecontrol = tune.control(), ...)
{
    call <- match.call()
    resp <- function(formula, data) model.response(model.frame(formula,
        data))
    classAgreement <- function(tab) {
        n <- sum(tab)
        if (!is.null(dimnames(tab))) {
            lev <- intersect(colnames(tab), rownames(tab))
            p0 <- sum(diag(tab[lev, lev]))/n
        }
        else {
            m <- min(dim(tab))
            p0 <- sum(diag(tab[1:m, 1:m]))/n
        }
        p0
    }
    if (!is.character(method))
        method2 <- deparse(substitute(method))
    else
        method2 <- method
    if (tunecontrol$sampling == "cross")
        validation.x <- validation.y <- NULL
    useFormula <- is.null(train.y)
    if (useFormula && (is.null(data) || length(data) == 0))
        data <- model.frame(train.x)
    if (is.vector(train.x))
        train.x <- t(t(train.x))
    if (is.data.frame(train.y))
        train.y <- as.matrix(train.y)
    if (!is.null(validation.x))
        tunecontrol$fix <- 1
    n <- nrow(if (useFormula)
        data
    else train.x)
    perm.ind <- sample(n)
    if (tunecontrol$sampling == "cross") {
        if (tunecontrol$cross > n)
            stop("`cross' must not exceed sampling size!")
        if (tunecontrol$cross == 1)
            stop("`cross' must be greater than 1!")
    }
    train.ind <- if (tunecontrol$sampling == "cross")
        tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
    else if (tunecontrol$sampling == "fix")
        list(perm.ind[1:trunc(n * tunecontrol$fix)])
    else lapply(1:tunecontrol$nboot, function(x) sample(n, n *
        tunecontrol$boot.size))
    parameters <- if (is.null(ranges))
        data.frame(dummyparameter = 0)
    else expand.grid(ranges)
    p <- nrow(parameters)
    if (!is.logical(tunecontrol$random)) {
        if (tunecontrol$random < 1)
            stop("random must be a strictly positive integer")
        if (tunecontrol$random > p)
            tunecontrol$random <- p
        parameters <- parameters[sample(1:p, tunecontrol$random),
            ]
    }
    model.errors <- c()
    for (para.set in 1:p) {
        sampling.errors <- c()
        for (sample in 1:length(train.ind)) {
            repeat.errors <- c()
            for (reps in 1:tunecontrol$nrepeat) {
                pars <- if (is.null(ranges))
                  NULL
                else lapply(parameters[para.set, , drop = FALSE],
                  unlist)
                model <- if (useFormula)
                  do.call(method, c(list(train.x, data = data,subset = train.ind[[sample]]), pars, list(...)))
                else do.call(method, c(list(train.x[train.ind[[sample]],], y = train.y[train.ind[[sample]]]), pars,list(...)))
                pred <- predict.func(model, if (!is.null(validation.x))
                  validation.x
                else if (useFormula)
                  data[-train.ind[[sample]], , drop = FALSE]
                else train.x[-train.ind[[sample]], , drop = FALSE])
                true.y <- if (!is.null(validation.y))
                  validation.y
                else if (useFormula)
                  resp(train.x, data[-train.ind[[sample]], ])
                else train.y[-train.ind[[sample]]]
                repeat.errors[reps] <- if (is.factor(true.y))
                  1 - classAgreement(table(pred, true.y))
                else crossprod(pred - true.y)/length(pred)
            }
            sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
        }
        model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
    }
    best <- which.min(model.errors)
    pars <- if (is.null(ranges))
        NULL
    else lapply(parameters[best, , drop = FALSE], unlist)
    structure(list(best.parameters = parameters[best, , drop = FALSE],
        best.performance = model.errors[best], method = method2,
        nparcomb = nrow(parameters), train.ind = train.ind, sampling = switch(tunecontrol$sampling,
            fix = "fixed training/validation set", bootstrap = "bootstrapping",
            cross = if (tunecontrol$cross == n) "leave-one-out" else paste(tunecontrol$cross,
                "-fold cross validation", sep = "")), performances = if (tunecontrol$performances) cbind(parameters,
            error = model.errors), best.model = if (tunecontrol$best.model) {
            modeltmp <- if (useFormula) do.call(method, c(list(train.x,
                data = data), pars, list(...))) else do.call(method,
                c(list(x = train.x, y = train.y), pars, list(...)))
            call[[1]] <- as.symbol("best.tune")
            modeltmp$call <- call
            modeltmp
        }), class = "tune")
}

  
  
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
 set.seed(rand)


 the.vector.of.all.parameters <- vector(length=0, mode="list")
 eset.matrix                  <- exprs(eset)  # row=genes, col=samples
 Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model






 train.matrix <- t(eset.matrix) #row=samples,col=genes
 train.factor <- class.column.factor


# now the 'best' parameter are calculated if this is nessessary ( if there are choices)
     if (! (all(sapply(poss.parameters,length) %in% 1))) {
    result <- Mytune(the.function.for.classification,
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
  
  predict.function <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))
     

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

ClassifierBuild.exprSetRG<- function(eset,
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

{ require(arrayMagic)
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
    warning("There is at least one parameter in the parameter list that does not correspond to a parameter of the classification, cluster or selection function.")

  if ("theParameter" %in% names(poss.parameters))
    stop("You have to change the parameter because 'theParameter','cluster.check' and 'throw.away' are not allowed as parameter names.")
  
  

  class.column.factor <- class.factor.format(getExprSetGreenMinusRed(eset),class.column,reference.class)
  levels.class        <- levels(class.column.factor)
  nlevels.class       <- nlevels(class.column.factor)
  nn                  <- length(class.column.factor)

  
  we                        <- new.env()
  #original.tune.environment <- environment(tune)
  #environment(tune)         <- we

  ## x is a matrix with rows=samples, columns=genes (for tune)
  ## m is a matrix with columns=samples, rows=genes (as usual with exprSets)
  tIp <- vector(mode="list",length=length(thePreprocessingMethods))




Mytune<-function (method, train.x, train.y = NULL, data = list(), validation.x = NULL,
    validation.y = NULL, ranges = NULL, predict.func = predict,
    tunecontrol = tune.control(), ...)
{
    call <- match.call()
    resp <- function(formula, data) model.response(model.frame(formula,
        data))
    classAgreement <- function(tab) {
        n <- sum(tab)
        if (!is.null(dimnames(tab))) {
            lev <- intersect(colnames(tab), rownames(tab))
            p0 <- sum(diag(tab[lev, lev]))/n
        }
        else {
            m <- min(dim(tab))
            p0 <- sum(diag(tab[1:m, 1:m]))/n
        }
        p0
    }
    if (!is.character(method))
        method2 <- deparse(substitute(method))
    else
        method2 <- method
    if (tunecontrol$sampling == "cross")
        validation.x <- validation.y <- NULL
    useFormula <- is.null(train.y)
    if (useFormula && (is.null(data) || length(data) == 0))
        data <- model.frame(train.x)
    if (is.vector(train.x))
        train.x <- t(t(train.x))
    if (is.data.frame(train.y))
        train.y <- as.matrix(train.y)
    if (!is.null(validation.x))
        tunecontrol$fix <- 1
    n <- nrow(if (useFormula)
        data
    else train.x)
    perm.ind <- sample(n)
    if (tunecontrol$sampling == "cross") {
        if (tunecontrol$cross > n)
            stop("`cross' must not exceed sampling size!")
        if (tunecontrol$cross == 1)
            stop("`cross' must be greater than 1!")
    }
    train.ind <- if (tunecontrol$sampling == "cross")
        tapply(1:n, cut(1:n, breaks = tunecontrol$cross), function(x) perm.ind[-x])
    else if (tunecontrol$sampling == "fix")
        list(perm.ind[1:trunc(n * tunecontrol$fix)])
    else lapply(1:tunecontrol$nboot, function(x) sample(n, n *
        tunecontrol$boot.size))
    parameters <- if (is.null(ranges))
        data.frame(dummyparameter = 0)
    else expand.grid(ranges)
    p <- nrow(parameters)
    if (!is.logical(tunecontrol$random)) {
        if (tunecontrol$random < 1)
            stop("random must be a strictly positive integer")
        if (tunecontrol$random > p)
            tunecontrol$random <- p
        parameters <- parameters[sample(1:p, tunecontrol$random),
            ]
    }
    model.errors <- c()
    for (para.set in 1:p) {
        sampling.errors <- c()
        for (sample in 1:length(train.ind)) {
            repeat.errors <- c()
            for (reps in 1:tunecontrol$nrepeat) {
                pars <- if (is.null(ranges))
                  NULL
                else lapply(parameters[para.set, , drop = FALSE],
                  unlist)
                model <- if (useFormula)
                  do.call(method, c(list(train.x, data = data,subset = train.ind[[sample]]), pars, list(...)))
                else do.call(method, c(list(train.x[train.ind[[sample]],], y = train.y[train.ind[[sample]]]), pars,list(...)))
                pred <- predict.func(model, if (!is.null(validation.x))
                  validation.x
                else if (useFormula)
                  data[-train.ind[[sample]], , drop = FALSE]
                else train.x[-train.ind[[sample]], , drop = FALSE])
                true.y <- if (!is.null(validation.y))
                  validation.y
                else if (useFormula)
                  resp(train.x, data[-train.ind[[sample]], ])
                else train.y[-train.ind[[sample]]]
                repeat.errors[reps] <- if (is.factor(true.y))
                  1 - classAgreement(table(pred, true.y))
                else crossprod(pred - true.y)/length(pred)
            }
            sampling.errors[sample] <- tunecontrol$repeat.aggregate(repeat.errors)
        }
        model.errors[para.set] <- tunecontrol$sampling.aggregate(sampling.errors)
    }
    best <- which.min(model.errors)
    pars <- if (is.null(ranges))
        NULL
    else lapply(parameters[best, , drop = FALSE], unlist)
    structure(list(best.parameters = parameters[best, , drop = FALSE],
        best.performance = model.errors[best], method = method2,
        nparcomb = nrow(parameters), train.ind = train.ind, sampling = switch(tunecontrol$sampling,
            fix = "fixed training/validation set", bootstrap = "bootstrapping",
            cross = if (tunecontrol$cross == n) "leave-one-out" else paste(tunecontrol$cross,
                "-fold cross validation", sep = "")), performances = if (tunecontrol$performances) cbind(parameters,
            error = model.errors), best.model = if (tunecontrol$best.model) {
            modeltmp <- if (useFormula) do.call(method, c(list(train.x,
                data = data), pars, list(...))) else do.call(method,
                c(list(x = train.x, y = train.y), pars, list(...)))
            call[[1]] <- as.symbol("best.tune")
            modeltmp$call <- call
            modeltmp
        }), class = "tune")
}


  

  preprocessing <- function(x,label,PreprocessingSlots=tIp,...)
   { m <- t(x)
     for (i in 1:length(thePreprocessingMethods))
       { gene.reduce             <- get(thePreprocessingMethods[i])(m,label,theParameter=PreprocessingSlots[[i]],,...)
         m                       <- gene.reduce$matrix
         PreprocessingSlots[[i]] <- gene.reduce$parameter
       }
    return(list(transf.matrix=t(m),ThePreSlots=PreprocessingSlots))
   }

  
  the.function.for.classification <-   function(x,y,...)
          {preprocessing.step <- preprocessing(x,y,...)
           tf.matrix          <- preprocessing.step$transf.matrix
           a                  <- ncol(tf.matrix)
           tf.matrix          <- tf.matrix[,1:(a/2)] - tf.matrix[,(a/2+1):a]  # now we have the log ratios for further things
           PreproParameter    <- preprocessing.step$ThePreSlots
           aaa                <- get(classification.fun)(tf.matrix,y,...)
           predict.function   <- aaa$predict
           assign("model.information",aaa$info,envir=we)

           rfct <- function(test.matrix)
           {testmatrix <- preprocessing(test.matrix,y,PreprocessingSlots=PreproParameter)$transf.matrix
            a          <- ncol(testmatrix)
            testmatrix <- testmatrix[,1:(a/2)] - testmatrix[,(a/2+1):a]
            result     <- predict.function(testmatrix)
            return(result)
           }
          return(rfct)
          }


  
 
  

  predicted <- function(model,test) return(model(test))


 ###### the random generator is set 
 set.seed(rand)



 the.vector.of.all.parameters <- vector(length=0, mode="list")
 eset.matrix              <- rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset)))  # row=genes, col=samples genes are two times
 Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model


 train.matrix <- t(eset.matrix) #row=samples,col=genes
 train.factor <- class.column.factor

  if (! (all(sapply(poss.parameters,length) %in% 1))){
     result<-Mytune(the.function.for.classification,
	          train.matrix,
	          train.factor,
	          ranges=poss.parameters,
	          predict.func=predicted,
	          tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
     parameter.list <- result$best.parameter
     the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
 } else {
     parameter.list <- as.data.frame(poss.parameters)
     }
 predict.function <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))
     

 classificationFunction <- function (eset){
   new.matr <- t(rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset))))
   if (ncol(train.matrix) != ncol(new.matr)) stop("The number of genes in the eset does not match the number of genes that are needed for the classifier.")
   return(predict.function(new.matr))
 }


 if (length(the.vector.of.all.parameters)!=0)
  {parameter.vector             <- unlist (the.vector.of.all.parameters)
   the.vector.of.all.parameters <- tapply(as.vector(parameter.vector),factor(names(parameter.vector)),c)
  }



## all data are collected
if(information)  
 Model.information.List[[1]] <- with(we,model.information)

  
  rv <- list(classifier=predict.function,
             classifier.for.exprSetRG= classificationFunction,
             classes=class.column.factor,
             parameter=the.vector.of.all.parameters,
             thePreprocessingMethods=thePreprocessingMethods,
             class.method=classification.fun,
             ross.inner=cross.inner,
             information=Model.information.List)

  #environment(tune)  <-  original.tune.environment

  return(rv)
}
