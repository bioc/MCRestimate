intersectList <- function (x)
  { if (class(x) != "list") stop("The first argument must be a list.")
  if (is.null(names(x)))
    gg <- paste("x[",1:length(x),"]",sep="")
    else{g <-which(names(x)=="")
         gg <- names(x)
         gg[g] <-  g
       }
    a <- length(x)
    b <- list()
    b[[1]] <- x[[1]]
    names(b) <- gg[1]
    for (i in 2:a){
      d <- lapply(b,function(y) intersect(y,x[[i]]))
      new.name <- gg[i]
      names(d) <- paste(names(d), "||", new.name)
      e <- list(unique(x[[i]]))
      names(e) <- new.name
      b <- c(d,b,e)
    }
  return(b)}
      
    
