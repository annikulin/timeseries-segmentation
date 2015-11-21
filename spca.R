spca <- function(a, b) {
  #Krzanowski-similarity between a and b hyperplanes
  trans_a <- t(a);
  trans_b <- t(b);
  
  m_trace <- sum(diag(trans_a %*% b %*% trans_b %*% a));
  ks <- m_trace / dim(a)[2];
  
  return(ks)
}
