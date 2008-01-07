###########################################################################################
##### Funktions for the vignette 'Basic compendium of computational diagnostic tools'######
#####                     Version 1.0                                                ######
#####  written by Markus Ruschhaupt and Ulrich Mannsmann  and Wolfgang Huber         ######
###########################################################################################


class.factor.format <-  function(x, class.column, reference.class=NULL)
 { #if (!"exprSet" %in% class(x))
    # stop ("'x' must be an exprSet")
   if (length(class.column)!=1)
     stop("'class.column' must have length 1.")
   if (is.character(class.column)) {
     if (!class.column %in% names(pData(x)))
       stop("The value of 'class.column' should be the name of a column of pData(x).")
   } else if (is.numeric(class.column)) {
     if ( class.column > ncol(pData(x)) | class.column < 0)
       stop("The value of 'class.column' should be between 1 and ncol(pData(x)).")
   } else {
     stop("'class.column' should be a character or numeric.")
   }
   class.vector <- pData(x)[,class.column]

   if (! (is.null(reference.class))) {       
     if (! all(reference.class %in% class.vector) ) 
       stop(paste("The values of 'reference.class' should be contained in the values of column", class.column, "of pData(x)."))
     if( length(reference.class) == 1 ){
       reference <- as.character(reference.class)
     }else{
       reference <- "Reference"
     }
     class.vector <- ifelse(class.vector %in% reference.class, 
                            reference, "Alternativ")
   }
   return(as.factor(class.vector))
 }
