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
    warning("There may be at least one parameter in the parameter list that does not correspond to a parameter of the classification or any preprocessing function.")

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
 

  the.function.for.classification<-  function(x,y,...)
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

 #######creating the block of samples which should be in the test set


 if (cross.outer > nn | cross.outer < 2 ) stop(paste("You have to specify at least two different groups and not more than",nn))
 if (cross.repeat < 1)        stop("There must be at least one repetition of the cross validation procedure")

 the.cut <- cut(1:nn,cross.outer,1:cross.outer)

 votal.matrix                 <- matrix(0,ncol=nlevels.class, nrow=nn)
 the.vector.of.all.parameters <- vector(length=0, mode="list")

 eset.matrix                  <- exprs(eset)  # row=genes, col=samples

  #assign("model.information",list(),envir=we)
  Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model



  for (l in 1:cross.repeat)
  { the.votes.per.cv             <- matrix(NA, ncol=nlevels.class,nrow=nn)
    rownames(the.votes.per.cv)   <- class.column.factor
    colnames(the.votes.per.cv)   <- levels.class
    if (stratify)
      sampleOfFolds  <- get("balanced.folds",en=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
    else
      permutated.cut <- sample(the.cut,nn)
    
  ## the start of the outer cross-validation ##
  #############################################

    for (sample in 1:cross.outer)
     {

   # dividing the set into a testset and a trainingsset.
    if (stratify){
      block <- rep(FALSE,nn)
      block[sampleOfFolds[[sample]]] <- TRUE
    }else{
     block <- permutated.cut==sample
    }
     train.matrix <- t(eset.matrix[,!block,drop=FALSE]) #row=samples,col=genes
     test.matrix  <- t(eset.matrix[,block,drop=FALSE])
     train.factor <- class.column.factor[!block]


   # now the 'best' parameter are calculated if this is nessessary ( if there are choices)
     if (! (all(sapply(poss.parameters,length) %in% 1)))
      {#if(length(train.factor)%/%cross.inner <= 1) stop(paste("Please choose a value for 'cross.inner' with the following attribute:",length(train.factor),"/ 'cross.inner' > 2"))
        result <- Mytune(the.function.for.classification,
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
  votal.matrix <- votal.matrix + the.votes.per.cv
  }

  votal.matrix <- votal.matrix /cross.repeat

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
 


## if plot.label is specified the vote matrix get new rownames
  
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
                                information=Model.information.List)
                           
   
  class(rv) <- "MCRestimate"
  #environment(tune)  <-  original.tune.environment
  
  return(rv)
}


## Version for exprSetRG (arrayMagic)

MCRestimate.exprSetRG<- function(eset,
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
    warning("There may be at least one parameter in the parameter list that does not correspond to a parameter of the classification, cluster or selection function.")

  if ("theParameter" %in% names(poss.parameters))
    stop("You have to change the parameter because 'theParameter' is not allowed as parameter names.")
  
  
  GreenMinusRedSet    <- getExprSetGreenMinusRed(eset) # it is not important if this is green minus red or red minus green because we are only interested in the phenoData.
  class.column.factor <- class.factor.format(GreenMinusRedSet,class.column,reference.class)
  levels.class        <- levels(class.column.factor)
  nlevels.class       <- nlevels(class.column.factor)
  nn                  <- length(class.column.factor)

  
  we                        <- new.env()
  #original.tune.environment <- environment(tune)
  #environment(tune)         <- we


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


  the.function.for.classification<-  function(x,y,...)
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
 
 #######creating the block of samples which should be in the test set
 

 if (cross.outer > nn | cross.outer < 2 ) stop(paste("You have to specify at least two different groups and not more than",nn))
 if (cross.repeat < 1)        stop("There must be at least one repetition of the cross validation procedure")

 the.cut <- cut(1:nn,cross.outer,1:cross.outer)

 votal.matrix                 <- matrix(0,ncol=nlevels.class, nrow=nn)
 the.vector.of.all.parameters <- vector(length=0, mode="list")

 eset.matrix              <- rbind(exprs(getExprSetGreen(eset)),exprs(getExprSetRed(eset)))  # row=genes, col=samples genes are two times

 #assign("model.information",list(),envir=we) 
 Model.information.List <- list() # a list that will contain additional model information like the genes that are used for the model



  
  for (l in 1:cross.repeat)
  { the.votes.per.cv             <- matrix(NA, ncol=nlevels.class,nrow=nn)
    rownames(the.votes.per.cv)   <- class.column.factor
    colnames(the.votes.per.cv)   <- levels.class
    if (stratify)
      sampleOfFolds  <- get("balanced.folds",en=asNamespace("pamr"))(class.column.factor, nfolds=cross.outer)
    else
      permutated.cut <- sample(the.cut,nn)
    
  ## the start of the outer cross-validation ##
  #############################################
    
    for (sample in 1:cross.outer)
     {

   # dividing the set into a testset and a trainingsset.
     if (stratify){
      block <- rep(FALSE,nn)
      block[sampleOfFolds[[sample]]] <- TRUE
    }else{
     block <- permutated.cut==sample
    }
     train.matrix <- t(eset.matrix[,!block,drop=FALSE]) #row=samples,col=genes
     test.matrix  <- t(eset.matrix[,block,drop=FALSE])
     train.factor <- class.column.factor[!block]


   # now the 'best' parameter are calculated if this is nessessary ( if there are choices)
     
     if (! (all(sapply(poss.parameters,length) %in% 1)))
      {#if(length(train.factor)%/%cross.inner <= 1) stop(paste("Please choose a value for 'cross.inner' with the following attribute:",length(train.factor),"/ 'cross.inner' > 2"))
        result <- Mytune(the.function.for.classification,
                         train.matrix,
		         train.factor,
                         ranges=poss.parameters,
                         predict.func=predicted,
                         tunecontrol=tune.control(cross=cross.inner,best.model=FALSE))
       parameter.list <- result$best.parameter
       the.vector.of.all.parameters <- c(the.vector.of.all.parameters,parameter.list)
       #predict.function <- result$best.model
     }else{
       parameter.list <- as.data.frame(poss.parameters)
     }
    predict.function         <- do.call("the.function.for.classification", c(list(x=train.matrix,y=train.factor),parameter.list))


   # the votes for a class are calculated
   if(information)
     Model.information.List[[(l-1)*cross.outer + sample]] <- with(we,model.information)
     
     pred.vector              <- predict.function(test.matrix)
     vote.matrix              <- t(sapply(1:length(pred.vector), function(j) as.numeric(levels.class==pred.vector[j])))
     colnames(vote.matrix)    <- levels.class
     the.votes.per.cv[block,] <- vote.matrix

   }
  votal.matrix <- votal.matrix + the.votes.per.cv
  }

  votal.matrix <- votal.matrix /cross.repeat

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
 


## if plot.label is specified the vote matrix get new rownames
  
  if (!(is.null(plot.label)))
    if (length(plot.label)==1 ) sample.names <- pData(GreenMinusRedSet)[,plot.label]
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
                                information=Model.information.List)
                           
   
  class(rv) <- "MCRestimate"
 # environment(tune)  <-  original.tune.environment
  
  return(rv)
}
