

MCRconfusion <- function(x, col.names=names(x),row.names=NULL)
  { if (class(x) != "list") stop("The first argument must be a list.")
    if (!(all(sapply(x,function(y) class(y))=="MCRestimate"))) stop("Each object in your given list must be a member of class 'MCRestimate'.")
    if (!is.null(col.names) & length(x) != length(col.names)) stop("The length of your list and the length of col.names must be equal.")
    ref.class <- x[[1]]$classes
    if (!(all(sapply(x,function(y) all(y$classes==ref.class))))) stop("At least one class lable does not fit.")
    number <- nlevels(ref.class)
    names <- levels(ref.class)
    reference <- table(ref.class)

    res.table <- matrix(NA, nrow=number+1,ncol=length(x)+1)
    if(is.null(col.names))
       colnames(res.table) <- c(paste("Method", 1:length(x)),"Size") else
       colnames(res.table) <- c(col.names,"Size")

    if (is.null(row.names))row.names <- names
    if (length(row.names)!=number) stop("There must be a row name for each sample group.")
    rownames(res.table) <- c(row.names,"All")

    for (i in 1:number)
     {x1 <- sapply(x, function(y) {a <- y$table ; sum(a[rownames(a) == names[i],!(colnames(a) %in% c(names[i],"class error"))])})
      x2 <- as.vector(reference[names(reference)==names[i]])
      res.table[i,] <-c(x1,x2)}
    res.table[number+1,] <- apply(res.table[1:number,],2,sum)                                                   
    return(res.table)}
