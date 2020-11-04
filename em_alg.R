

#' Calculates mij for the E step
#'
#' @param proj list of information for this projection
#' @param theta matrix of estimated theta values
#' @param j pixel of interest
#' @return mij
#' @export
mij <- function(proj, theta, j) {
  
  d <- proj$d
  
  idx <- proj$idx
  
  y <- proj$y
  
  
  
  if (j %in% proj$idx) {
    
    k <- which(proj$idx == j)
    
  } else {
    
    return(0)
    
  }
  
  
  
  if (k > 1) {
    
    idx_sub <- idx[1:(k - 1)]
    
    theta_sub <- theta[idx_sub]
    
  } else {
    
    theta_sub <- 0
    
  }
  
  return(d * (exp(-sum(theta_sub)) - exp(-sum(theta[idx]))) + y)
  
}
  
#' Calculates nij for the E step
#'
#' @param proj list of information for this projection
#' @param theta matrix of estimated theta values
#' @param j pixel of interest
#' @return nij
#' @export
nij <- function(proj, theta, j) {
  
  d <- proj$d
  
  idx <- proj$idx
  
  y <- proj$y
  
  
  
  if (j %in% proj$idx) {
    
    k <- which(proj$idx == j)
    
  } else {
    
    return(0)
    
  }
  
  
  
  idx_sub <- idx[1:k]
  
  theta_sub <- theta[idx_sub]
  
  
  
  return(d * (exp(-sum(theta_sub)) - exp(-sum(theta[idx]))) + y)
  
}
  



  
  
  