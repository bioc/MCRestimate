print.MCRestimate <- function(x,...)
 { class.function        <- x$class.method
   type.of.cv <- ifelse(x$cross.outer==length(x$classes),"leave-one-out",paste(x$cross.outer,"-fold",sep=""))
   one.or.many <- ifelse(x$cross.repeat==1,paste(x$cross.repeat,"repetition of"),paste(x$cross.repeat,"repetitions of"))
   thePres <- x$thePreprocessingMethods
   thePres <- thePres[thePres!="identity"]
   nP <- length(thePres)
   cat("\n")
   cat("Result of MCRestimate with", one.or.many,type.of.cv,"cross-validation",fill=TRUE)
   cat("\n")
   if(nP>0){
     for (i in 1:nP) cat(paste("Preprocessing function",i,":",thePres[i]),"\n")
   } else
     cat("No preprocessing \n")
     cat( "Classification function  : ", x$class.method, "\n\n",
       "The confusion table: \n",sep="",...)
   print(round(x$table,3), ...)
 }
