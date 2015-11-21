pcaresid <- function(x,ndim,flag) {
  #Computes the cost of a segment
  #Inputs: 
  #         x: time series to be segmented (columns: variables)
  #         ndim: number of retained principal components
  #         flag: cost type (0: Hotelling (T2), 1: avarage residual error(Q))
  # Outputs
  #         pc: retained principal components of the segment
  #         cost: cost of the segment
  #         avg: mean of the segment

  m <- dim(x)[1];
  n <- dim(x)[2];
  
  if (is.null(n)) {
    n <- 1;
  }
  
  if (length(dim(ndim)) > 1) {
    stop('ndim must be a scalar value!');
  }
  
  if (ndim >= n) {
    stop('ndim must be smaller than the variables (columns) of x!');    
  }

  pca <- prcomp(x); # execute principal component analysis
  avg <- colMeans(x); # get the average of x for each variable (columns)
  avgx <- matrix(rep(avg, m), nrow = m, byrow = T); # just repeat the avarage for each variable for further computation
  retain <- t(pca$rotation[,1:ndim]); # get the retained principal components and transfer their matrices for further usage
  predictx <- avgx + pca$x[,1:ndim] %*% retain; # predicted values of x using the retained principal components
  residuals <- x - predictx; # get the prediction error for each time stamp
  
  if (flag == 1) {
    cost <- mean(rowSums(residuals ^ 2)); # average residual error(Q)    
  } else {
    cost <- 0; # Hotelling (T2)
  }
  
  pc <- pca$rotation[,1:ndim]; # retained principal components of the segment
  
  returnlist <- list(cost=cost, pc=pc, avg=avg);
  return(returnlist);
  #return(cost);
}
