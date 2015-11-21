dtwpca <- function(seqx,seqy) {
  # Dynamic time warping of previously segmented multivariate time series using 
  # the retaind principal components of each segment. The basic distance
  # between the segments is 1 - Kranowski similairty of the retained
  # principal components
  # Note: DTW is not optimized
  # Inputs:
  #         seqx, seqy: arrays of structure, describe the segmented time series
  # Output:     
  #         dist = DTW distance of the two multivariate time series

  N <- length(seqx);
  M <- length(seqy);
  
  # initialize the cummulated distance matrix
  D <- matrix(data = rep(Inf, times = M*N), nrow = M, ncol = N);
  
  # compute dtw
  D[1,1] <- 1 - spca(seqx[[1]]$pc, seqy[[1]]$pc); # creates distance from Krzanowski-similarity
  for (i in 2 : M) {
    D[i,1] <- D[i-1,1] + 1 - spca(seqy[[i]]$pc, seqx[[1]]$pc);    
  }

  for (i in 2 : N) {
    D[1,i] <- D[1,i-1] + 1 - spca(seqy[[1]]$pc, seqx[[i]]$pc);
  }
  
  
  for (n in 2 : N) {
    for (m in 2 : M) {
      D[m,n] = min(D[m, n-1], D[m-1, n], D[m-1, n-1]) + 1 - spca(seqx[[n]]$pc,seqy[[m]]$pc);
    }
  }
  
  dist = D[M,N];
  return(dist);
}