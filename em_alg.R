

#' Calculates mij for the E step
#'
#' @param proj list of information for this projection
#' @param theta matrix of initial theta values
#' @param j pixel of interest
#' @return mij
#' @export
mij <- function(proj, theta, j) {
  
  d <- proj$d
  
  idx <- proj$idx
  
  y <- proj$y
  
  l <- proj$l
  
  
  
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
  
  return(d * (exp(-sum(l * theta_sub)) - exp(-sum(l * theta[idx]))) + y)
  
}



#' Calculates nij for the E step
#'
#' @param proj list of information for this projection
#' @param theta matrix of initial theta values
#' @param j pixel of interest
#' @return nij
#' @export
nij <- function(proj, theta, j) {
  
  d <- proj$d
  
  idx <- proj$idx
  
  y <- proj$y
  
  l <- proj$l
  
  
  
  if (j %in% proj$idx) {
    
    k <- which(proj$idx == j)
    
  } else {
    
    return(0)
    
  }
  
  
  
  idx_sub <- idx[1:k]
  
  theta_sub <- theta[idx_sub]
  
  
  
  return(d * (exp(-sum(l * theta_sub)) - exp(-sum(l * theta[idx]))) + y)
  
}



#' Maximizes Q function for pixel j to estimate theta j
#'
#' @param theta_j theta for pixel j to optimize for
#' @param proj_list list of projections
#' @param theta matrix of initial theta values
#' @param j pixel of interest
#' @return q function value for theta j
#' @export
q_fun_j <- function(thetaj, proj_list, theta, j) {
  
  num_proj <- length(proj_list)
  
  val <- 0
  
  for (i in 1:num_proj) {
    
    l <- proj_list[[i]]$l
    
    m_exp <- mij(proj_list[[i]], theta, j)
    
    n_exp <- nij(proj_list[[i]], theta, j)
    
    val <- val - n_exp * l + (m_exp - n_exp) * l / (exp(l * thetaj) - 1)
    
  }
  
  return(val)
  
}



#' Detects whether there are any projections with y = 0 in proj_list
#'
#' @param proj_list list of projections
#' @return TRUE/FALSE depending if a zero was detected or not
#' @export
y_zero <- function(proj_list) {
  
  has_zero <- FALSE   
  
  for (proj in proj_list) {
    
    if (proj$y == 0) {
      has_zero <- TRUE
    }
    
  }
  
  return(has_zero)
  
}



# To bypass the theta = zero error
max_q_fun_j <- function(interval, proj_list, theta, j) {
  
  e <- try(
    # if I need more efficiency, calculate nij and mij outside of function
    theta_est <- uniroot(q_fun_j, interval, proj_list, theta, j, extendInt = "downX")$root, 
    silent = TRUE
  )
  
  
  if (class(e) == "try-error") {
    return(0)
  } else {
    return(theta_est)
  }
  
}



#' EM algorithm for transmission tomography
#'
#' @param proj_list list of projections
#' @param theta matrix of initial theta values. Assumes the -1's in the matrix refer to dead space to 
#' be ignored
#' @param tol tolerance used for the stopping rule
#' @return estimated theta values
#' @export
em_alg <- function(proj_list, theta, tol) {
  
  if(y_zero(proj_list)){
    return(simpleError("There's a projection with zero photons observed."))
  }
  
  theta_est <- matrix(-1, nrow = nrow(theta), ncol = ncol(theta))
  
  # indices of theta_est that are nonnegative
  nonneg_idx <- which(theta >= 0)
  
  ctr <- 0
  
  diff <- Inf
  
  while (diff > tol & ctr <= 100000) {
    
    for (j in nonneg_idx) {
      
      theta_est[j] <- max_q_fun_j(interval = c(0.0000001, 10), proj_list, theta, j)
      
    }
    
    diff <- sum((theta_est - theta)^2)
    
    theta <- theta_est
    
    ctr <- ctr + 1
    
  }
  
  return(list(theta_est = theta_est, ctr = ctr))
  
}






