pcaseg <- function(data,num_segments,q,flag) {
  # Bottom-up segmentation of multivariate time series using PCA based cost functions 
  # Inputs:
  #    data: multivariate time series (columns: variables)
  #    num_segments: desired number of segments
  #    q: number of principal compnents to be retained
  #    flag: cost type (0: Hotelling (T2), 1: avarage residual error(Q))
  
  print(Sys.time());
  
  minres <- ceiling(dim(data)[1]/600); # initial segmentation -> MUST BE MODIFIED BASED ON YOUR DATA!!!
  left_x <- seq(from = 1, to = dim(data)[1] - 1, by = minres); # starting points of the segments
  right_x <- left_x + minres; # ending points of the segments
  right_x[length(right_x)] <- dim(data)[1]; # last element is obviously is the size of the time series
  number_of_segments <- length(left_x); # initial number of segment

  # intialize segments
  segments <- list();
  for (i in 1 : number_of_segments) {
    segment <- NULL;
    segment$lx <- left_x[i];
    segment$rx <- right_x[i];
    segment$mc <- Inf;
    segment$c <- Inf;
    segments[[i]] <- segment;
  }
  
  tc <- list();
  
  # compute merge costs (i.e. cost of two consecutive segments)
  
  for (i in 1 : (number_of_segments - 1)) {
    sx <- data[segments[[i]]$lx : segments[[i+1]]$rx, ];
    segments[[i]]$mc = pcaresid(sx,q,flag)$cost; # compute merge cost with the consecutive segment
    sx <- data[segments[[i]]$lx : segments[[i]]$rx, ];
    segments[[i]]$c = pcaresid(sx,q,flag)$cost; # cost of the segment itself
  }
  
  # special handling of the last segment -> no merge cost
  sx <- data[segments[[i+1]]$lx :segments[[i+1]]$rx, ];
  segments[[i+1]]$c <- pcaresid(sx,q,flag)$cost;

  # segments are merged until the desired segment number is not reached
  while (length(segments) > num_segments) {
    i <- which.min(lapply(segments, function(x) x$mc)); # get the minimum of the merge costs
    temp <- segments[[i]]$mc;
    
    if (i > 1 && i < length(segments) - 1) { 
      # special case 1: neither the first nor last segment is merged
      segments[[i]]$c <- segments[[i]]$mc; # the merge cost of the two segments now become the cost of the merged segment
      sx <- data[segments[[i]]$lx : segments[[i+2]]$rx, ];
      segments[[i]]$mc <- pcaresid(sx,q,flag)$cost; # compute the new merge cost of the currently merged segment and the next one
      segments[[i]]$rx <- segments[[i+1]]$rx; # update the last data point of the newly merged segment
      segments[[i+1]] <- NULL; # delete the second segments of two we just merged
      i <- i - 1; # decrease index
      sx <- data[segments[[i]]$lx : segments[[i+1]]$rx, ]; 
      segments[[i]]$mc <- pcaresid(sx,q,flag)$cost; # update the merge cost of the segment in front of the newly merged segment
    } else if (i == 1) {  
      # special case 2: first segment is merged
      segments[[i]]$c <- segments[[i]]$mc;
      sx <- data[segments[[i]]$lx :segments[[i+2]]$rx, ];
      segments[[i]]$mc <- pcaresid(sx,q,flag)$cost;
      segments[[i]]$rx <- segments[[i+1]]$rx;
      segments[[i+1]] <- NULL;
      sx <- data[segments[[i]]$lx : segments[[i]]$rx, ];
      segments[[i]]$c <- pcaresid(sx,q,flag)$cost;
    } else {
      # special case 3: last segment is merged
      segments[[i]]$rx <- segments[[i+1]]$rx;
      segments[[i]]$c <- segments[[i]]$mc;
      segments[[i]]$mc <- Inf;
      segments[[i+1]] <- NULL;
      i <- i - 1;
      sx <- data[segments[[i]]$lx : segments[[i+1]]$rx, ];
      segments[[i]]$mc <- pcaresid(sx,q,flag)$cost;
    }
    
    #tc <- Reduce("+", lapply(segment, function(x) x$c)); #[tc; sum([segment.c])];
  }

  # go trough all segments and compute their merge costs (won't be used
  # anymore), principal components (will be used to compute the distance of
  # the segments) and avarage of the segments (won't be used anymore)
  
  for (i in 1 : length(segments)) {
    sx <- data[segments[[i]]$lx : segments[[i]]$rx, ];
    pcares <- pcaresid(sx, q, flag);
    segments[[i]]$mc <- pcares$cost;
    segments[[i]]$pc <- pcares$pc;
    segments[[i]]$avg <- pcares$avg;
  }
  print(Sys.time());
  return(segments);
}



