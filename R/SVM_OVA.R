#SVM with 'One Versus All' multiclass approach
SVM.OVA.wrap <-
function (x, y, gamma = NULL, kernel = "radial", ...)
{
    require(e1071)
    nclasses <- length(levels(y))
    #model.list will contain one model per class; every model is a binary classifier distinguishing elements of one class versus all others
    model.list<- list()
    y.levels<-levels(y)
    for (m in 1:nclasses)  {
      a.y <- rep("c1",length(y))
      #select all elements in y of class levels(y)[m]
      a.y[as.character(y)==(y.levels[m])]<-"c2"
      if (!is.null(gamma)) {
        gamma <- gamma/ncol(x)
        a.model <- svm(x, factor(a.y), gamma = gamma, kernel = as.character(kernel), ...)
    }
    else a.model <- svm(x, factor(a.y), kernel = as.character(kernel), ...)
    
    model.list<-c(model.list, list(a.model))
    
    }
    #
    predict.function <- function(testmatrix) {
      p.m=matrix(0,nrow=dim(testmatrix)[1],ncol=nclasses)
      for (m in 1:nclasses)  {
        p=attr(predict(model.list[[m]],testmatrix,decision.values=T), "decision.values")
	#elements of the class which appears first in y get positive decision values
	#therefore, the decision values for the classifier which distinguishes elements 
	#of the class which appears first in y versus all other elements has to be multipied by -1
	if (y.levels[m]==as.character(y[1]))
	  p.m[,m]= -1*p
	else
	  p.m[,m]= p
      }
      res.p<-NULL
      for (i in 1:(dim(p.m)[1]))  {
        a.v<-p.m[i,]
	res.p<-c(res.p, y.levels[which.min(a.v)])
      }
      return(as.factor(res.p))
    }
    return(list(predict = predict.function, info = list()))
}
