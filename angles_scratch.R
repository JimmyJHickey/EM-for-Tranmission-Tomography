
# detecting if any projections have length sqrt(2) (or not 1), indicating non-90 degree angles

d_sqrt2 <- function(proj_list) {
  
  has_sqrt2 <- FALSE   
  
  for (proj in proj_list) {
    
    if (proj$l != 1) {
      has_sqrt2 <- TRUE
    }
    
  }
  
  return(has_sqrt2)
  
}

# looks like no 45 degree angles are making it



neg_idx <- function(proj_list) {
  
  has_neg <- FALSE   
  
  for (proj in proj_list) {
    
    if (length(which(proj$idx <= 0)) >= 1) {
      has_neg <- TRUE
    }
    
  }
  
  return(has_neg)
  
}

# I was worried there might be negative indices, but it's all good


