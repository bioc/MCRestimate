
###############################################################################
# the prediction function only works with two different classes             ###
###############################################################################


PLR<- function (trainmatrix, resultvector, kappa = 0, eps = 1e-04)
{   if (!length(resultvector) == nrow(trainmatrix))
        stop("The data matrix must have as many rows as the class label vector.")
    if (!length(unique(resultvector)) == 2)
        stop("PLR only works with exactly two different classes")
    class.factor <- as.factor(resultvector)
    factorlevel <- list(level0 = c(), level1 = c())
    factorlevel$level0 <- levels(class.factor)[1]
    factorlevel$level1 <- levels(class.factor)[2]
    resultvector <- as.numeric(class.factor) - 1
    m <- nrow(trainmatrix)
    n <- ncol(trainmatrix)
    trainmatrix.c <- apply(trainmatrix, 2, function(k) k - mean(k))
    s <- svd(trainmatrix.c, nu = m)
    R <- s$u * matrix(s$d, m, m, byrow = TRUE)
    ty <- rep(1, m)
    aics <- trs <- double()
    results <- list()
    for (i in 1:length(kappa)){
      kapp <- kappa[i]
      FIT     <- plr.fit(resultvector, ty, R, kapp, eps = eps)
      results[[i]] <- FIT
      print(c(log10(kapp), FIT$aic, FIT$tr))
    }
    aics <- sapply(results,function(x) x$aic)
    trs <- sapply(results,function(x) x$tr)
    aicbest <- which.min(aics)
    res <- results[[aicbest]]
    b <- s$v %*% res$a[-1]
    result <- list(a = res$a[1], b = b, factorlevel = factorlevel, aics=aics,trs=trs,kappas=kappa)
    class(result) <- "PLR"
    return(result)
}


predict.PLR <- function(object,...)
  { testmatrix <- list(...)[[1]]
    X1c        <- apply(testmatrix, 2, function(x) x-mean(x))
    eta1       <- object$a + X1c %*% object$b
    pp         <- 1/(1+exp(-eta1))
    result     <- ifelse(pp>0.5,object$factorlevel$level1,object$factorlevel$level0)
   return(result)
   }


######################################################
#this functions are needed for the function above   ##
######################################################

 plr.fit <- function(y, ty, X0, kappa=0, eps=1e-4)  
 {
  m <- nrow(X0)
  X0 <- apply(X0, 2, function(x) x-mean(x))
  X <- cbind(rep(1,m), X0)
  n <- ncol(X)
 
 #
  p0 <- (sum(y)+1)/(sum(ty)+2)
  a <- rep(0, n)
  a[1] <- log(p0/(1-p0))
  R <- diag(x=kappa, nrow=n)
  R[1,1] <- 0

 #
 repeatit <- 1
 it <- 0
 while(repeatit) {
    it <- it + 1
    z <- X %*% a
    a1 <- a
    p <- 1/(1+exp(-z))
    mu <- ty*p
    w <- mu * (1-p)
    W <- matrix(w, m, n)
    G0 <- t(X) %*% (W*X)
    G <- G0 + R
    a <- solve(G, t(X) %*% (y-mu+w*z))
    da <- max(abs(a-a1))
    repeatit <- (da / max(abs(a))) > eps
 }
 dev <- bindev(y, ty, p)
 tr <- sum(diag(solve(G)%*%G0))
 aic <- dev + 2 * tr
 return(list(p=p, a=a, aic=aic, tr=tr, z=z))
 }


 bindev <- function(y, ty, p)
 {
   mu <- ty * p
   r <- ty - y
   f1 <- (y > 0)
   f2 <- (y < ty)
   d1 <- sum(y[f1] * log(y[f1] / mu[f1]))
   d2 <- sum(r[f2] * log(r[f2] / (ty[f2] - mu[f2])))
   d = 2 * d1 + 2 * d2;
 return(d)
 }




