



################################################################################################################
########        for each new classification method you have to write a function ################################
# the arguments : x = trainmatrix (row=sample, col =gene)                                                      #
#                 y = factor-the length of the factor must be equal to the number of rows of the trainmatrix   #
#                 additional argments for the classification method                                            #             
#                                                                                                              #
# the result    : prediction function for the trained model                                                    #
################################################################################################################
################################################################################################################
  


###########################
##### RandomForest     ####
###########################


RF.wrap <- function (x,y,...)
  { require(randomForest)
    forest <- randomForest(x,y,importance=TRUE,...)
    names <- forest$importance[(forest$importance[,2] != 0),]
    names <- names[order(names[,2],decreasing=TRUE),]
    names <- cbind(rownames(names),names)
    predict.function <- function(testmatrix) return(predict(forest,testmatrix))

    return(list(predict=predict.function,info=names))
    }



###################################################
######### Penaliesed logistic regression ##########
###################################################


PLR.wrap <- function(x,y,kappa=0,eps=1e-4,...)
  { model <- PLR(x,y,kappa,eps)
    predict.function <- function(testmatrix) return(predict(model,testmatrix))
   
   return(list(predict=predict.function,info=list()))
 }


###################################
########## PAM ####################
###################################


PAM.wrap <- function(x,y,threshold,...)
  { require(pamr)
    if (is.null(colnames(x))) colnames(x) <- rep(1,ncol(x))
    trainlist<- list(x=t(x),y=y,geneid=colnames(x))
    model <- pamr.train(trainlist,threshold=threshold)
    if (model$nonzero != 0){
      ## we include the function form the pamr package here, but remove a bug which appears when there is only one gene left
      ## we do print the results on screen
      pamr.listgenes2 <- function(fit, data, threshold, genenames = FALSE){
       if (is.null(fit$newy)) {
        y <- factor(data$y[fit$sample.subset])
       }
       if (!is.null(fit$newy)) {
        y <- factor(fit$newy[fit$sample.subset])
       }
       x <- data$x[fit$gene.subset, fit$sample.subset]
       if (genenames) {
        gnames <- data$genenames[fit$gene.subset]
       }
       if (!genenames) {
        gnames <- NULL
       }
       geneid <- data$geneid[fit$gene.subset]
       nc <- length(unique(y))
       aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
       cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
       d <- (cen - fit$centroid.overall)[aa, ,drop=FALSE]/fit$sd[aa]
       oo <- order(-apply(abs(d), 1, max))
       d <- round(d, 4)
       g <- gnames[aa]
       g1 <- geneid[aa]
       if (is.null(gnames)) {
        gnhdr <- NULL
       }
       if (!is.null(gnames)) {
        gnhdr <- "name"
       }
       options(width = 500)
       schdr <- paste(dimnames(table(y))$y, "score", sep = " ")
       res <- cbind(as.character(g1), g, d)[oo, ,drop=FALSE] #changed
       dimnames(res) <- list(NULL, c("id", gnhdr, schdr))
       return(res) # change
       }
    names <- pamr.listgenes2(model,trainlist,threshold)
    }
    else names <- c()
    
    predict.function <- function(testmatrix) return(pamr.predict(model,t(testmatrix),threshold=threshold))      
    return(list(predict=predict.function,info=names))
  }



############################
########  SVM   ############
############################


SVM.wrap <- function(x,y,gamma=NULL,kernel="radial",...)
  { require(e1071)
    if (!is.null(gamma)) {gamma <- gamma/ncol(x);model <- svm(x,y,gamma=gamma,kernel=as.character(kernel),...)}
    else model <- svm(x,y,kernel=as.character(kernel),...)
    predict.function <- function(testmatrix) return(predict(model,testmatrix))
    
    return(list(predict=predict.function,info=matrix(colnames(x),ncol=1)))
  }




#######################################
## wrapper for GPLS (for two groups) ##
#######################################


GPLS.wrap <- function (x,y,...){
  library(gpls)
    level.y <- levels(y)
    if (length(level.y)!=2) stop("Up to now this methods only works with two groups")
    y <- as.integer(y) - 1 # because if a factor is converted into an integer it starts with 1
    model <- glpls1a(x,y)
    #if (!model$convergence) stop("One of the models did not converged! You may have to increase the number of interations")
    cat(".")
    predict.function <- function(testmatrix)
     { beta <- model$coefficients
       res <-getFromNamespace("glpls1a.predict","gpls")(testmatrix,beta)
       result <- factor(ifelse(res > 0.5,level.y[2],level.y[1]))
       return(result)
      }

    return(list(predict=predict.function,info=list(cv=model$convergence)))
    }


