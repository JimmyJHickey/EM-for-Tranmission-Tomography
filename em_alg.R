

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
  


#' Maximizes Q function for pixel j to estimate theta j
#'
#' @param theta_j theta for pixel j to optimize for
#' @param proj_list list of projections
#' @param theta matrix of estimated theta values
#' @param j pixel of interest
#' @return q function value for theta j
#' @export
q_fun_j <- function(thetaj, proj_list, theta, j) {
  
  num_proj <- length(proj_list)
  
  val <- 0
  
  for (i in 1:num_proj) {
    
    m_exp <- mij(proj_list[[i]], theta, j)
    
    n_exp <- nij(proj_list[[i]], theta, j)
    
    val <- val - n_exp + (m_exp - n_exp) / (exp(thetaj) - 1)
    
  }
  
  return(val)
  
}
  
  

#' EM algorithm for transmission tomography
#'
#' @param proj_list list of projections
#' @param theta matrix of initial theta values
#' @param tol tolerance used for the stopping rule
#' @return estimated theta values
#' @export
em_alg <- function(proj_list, theta, tol) {
  
  theta_est <- matrix(NA, nrow = nrow(theta), ncol = ncol(theta))
  
  num_pixel <- nrow(theta) * ncol(theta)
  
  ctr <- 0
  
  diff <- Inf
  
  while (diff > tol & ctr <= 1000) {
    
    for (j in 1:num_pixel) {
      
      # what are the bounds of theta?
      # if I need more efficiency, calculate nij and mij outside of function
      theta_est[j] <- uniroot(q_fun_j, interval = c(0, 10), proj_list, theta, j)$root  
      
    }
    
    diff <- sum((theta_est - theta)^2)
    
    theta <- theta_est
    
    ctr <- ctr + 1
    
  }
  
  return(list(theta_est = theta_est, ctr = ctr))
  
}





  