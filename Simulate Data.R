#ST793 Final Project
#EM algorithm for CT Data
#Eric Yanchenko (Alvin Sheng and Jimmy Hickey) [Copyright]
#November 4, 2020

# calculates column-majorized vector index of a matrix index
calc_index = function(nrow, x, y){
  return(nrow * (y-1) + x)
}



# inverse of calc_index (for testing purposes)
calc_matrix_index <- function(idx, num_row) {
  
  row_idx <- rep(NA, length(idx))
  
  col_idx <- rep(NA, length(idx))
  
  for (i in 1:length(idx)) {
    
    col_idx[i] <- ceiling(idx[i] / num_row) 
    
    row_idx[i] <- idx[i] - (col_idx[i] - 1) * num_row
    
  }
  
  return(cbind(row_idx, col_idx))
  
} 



# returns a logical vector indicating the entries that are -1 and at either end of the vector
# I use this method of taking out -1's, so I can detect weird cases in which 
# -1 is in the middle of the vector instead
neg_ends <- function(hit_thetas) {
  
  p <- length(hit_thetas)
  
  nends <- rep(FALSE, length(hit_thetas))
  
  i <- 1
  
  while (hit_thetas[i] < 0 && i <= p) {
    nends[i] <- TRUE
    i <- i + 1
  }
  
  if (i == p) { # if it went through entire vector already
    return(nends)
  }
  
  # go from the other side
  
  other_i <- p
  
  while(hit_thetas[other_i] < 0 && nends[other_i] == FALSE && other_i >= 1) {
    nends[other_i] <- TRUE
    other_i <- other_i - 1
  }
  
  return(nends)
  
}



#Function which generates the observed data for a single theta vector
data_gen <- function(theta, d, l){
  #Theta is the given theta vector used to compute attenutation probs.
  #d is mean of initial Poisson distribution
  #l is the length of each projection through each pixel
  
  p = length(theta)
  
  if(p<=0){#If empty theta vector is given, return NA
    return(NA)
  }
  
  X1 = rpois(1, lambda = d) #Initial number of counts from Poisson dist.
  
  X = c(X1, rep(0, p)) #Vector to hold the number of counts after each interaction
  
  for(j in 1:p){
    X[j+1] <- rbinom(1, size = X[j], prob = exp(-theta[j] * l))
    #Generate the number of photons from Binom(X_i, exp(-theta_i * l))
  }
  
  Y = X[p+1] #Actually observed data (number of counts)
  
  return(Y)
}


#' Function which generates the observed data for entire theta matrix and arbitrary projections
#'
#' @param THETA input theta matrix used to compute probs.
#' @param d Poisson mean for initial Poisson generation
#' @param ROW rows which we project on (vector) Enter negative row number to get opposite direction
#' @param COL columns which we project on (vector) Enter negative col number to get opposite direction
#' @param reps number of times to run the function
#' @param rise_vec vector of numerator(s) of the slopes of the parallel projections (projections perpendicular to them will 
#' also be generated)
#' @param run_vec vector of denominator(s) of the slopes of the parallel projections (projections perpendicular to them will 
#' also be generated)
#' Note: One of rise or run is assumed to be one, and the other is assumed to be an integer >= one.
#' @return estimated theta values
#' @export
data_gen_df <- function(THETA, d, ROW, COL, reps=1, rise_vec = 1, run_vec = 1){
  
  r = dim(THETA)[1] #Number of rows in THETA
  c = dim(THETA)[2] #Number of cols in THETA
  
  #Check that ROW and COL projection indices are valid
  if(max(abs(ROW))> r || max(abs(COL)) > c){
    return(simpleError("Row and/or column indices out of range."))
  }
  
  #Return projection list
  #with d, length of projection within each pixel, indices (in order) beam went through, counts
  proj.list <- vector("list", 1)
  
  pl_idx <- 1
  
  for(a in 1:reps){
    
    #Run over ROWS first
    for(rr in ROW){
      if(rr > 0){
        y <- data_gen(THETA[rr,THETA[rr,]>=0], d, l = 1)
        #Matrix indices that this beam goes through
        idx = seq(rr, rr+(c-1)*r, r)
        idx = idx[which(THETA[rr,]>=0, arr.ind=TRUE)] #Drop indices with negative thetas
      }else{ #switch order of theta if row number is negative
        y <- data_gen(rev(THETA[-(rr),THETA[-(rr),]>=0]), d, l = 1)
        #Matrix indices that this beam goes through
        idx = seq(-rr, -rr+(c-1)*r, r)
        idx = rev(idx[which(THETA[-rr,]>=0, arr.ind=TRUE)])
      }
      
      if (!is.na(y)) {
        add.list <- list(d = d, l = 1, idx = idx, y = y)
        proj.list[[pl_idx]] <- add.list
        pl_idx <- pl_idx + 1
      }
      
    }
    
    #Run over COLS second
    for(cc in COL){
      if(cc >0){
        y<- data_gen(THETA[THETA[,cc]>=0,cc], d, l = 1)
        idx = seq((cc-1)*r+1, cc*r, 1)
        idx = idx[which(THETA[,cc]>=0, arr.ind = TRUE)] # Drop indices with negative thetas
      }else{
        y <- data_gen(rev(THETA[THETA[,-(cc)]>=0,-(cc)]), d, l = 1)
        idx = seq((-cc-1)*r+1, -cc*r, 1)
        idx = rev(idx[which(THETA[,-cc]>=0, arr.ind=TRUE)])
      }
      
      if (!is.na(y)) {
        add.list <- list(d = d, l = 1, idx = idx, y = y)
        proj.list[[pl_idx]] <- add.list
        pl_idx <- pl_idx + 1
      }
      
    }
    
    #Check that rise_list and run_list have the same lengths
    if(length(rise_vec) != length(run_vec)){
      return(simpleError("rise_vec and run_vec have different numbers of elements"))
    }
    
    if (length(rise_vec) != 0) {
      for (i in 1:length(rise_vec)) {
        # call a separate function for the angular projections
        angular_proj_list <- angular_proj_list_gen(THETA, d, ROW, COL, rise = rise_vec[i], run = run_vec[i])
        proj.list <- c(proj.list, angular_proj_list)
        pl_idx + length(angular_proj_list)
      }
    }
    
  }
  
  #Return list with d, l, idx, y for each 
  return(proj.list)
}



#' Function to generate one set of angular projections, to be used within data_gen_df
#'
#' @param THETA input theta matrix used to compute probs.
#' @param d Poisson mean for initial Poisson generation
#' @param ROW rows which we project on (vector) Enter negative row number to get opposite direction
#' @param COL columns which we project on (vector) Enter negative col number to get opposite direction
#' @param rise numerator of the slope of the parallel projections (projections perpendicular to them will 
#' also be generated)
#' @param run denominator of the slope of the parallel projections (projections perpendicular to them will 
#' also be generated)
#' Note: One of rise or run is assumed to be one, and the other is assumed to be an integer >= one.
#' @return estimated theta values
#' @export
angular_proj_list_gen <- function(THETA, d, ROW, COL, rise, run) {
  
  r = dim(THETA)[1] #Number of rows in THETA
  c = dim(THETA)[2] #Number of cols in THETA
  
  rise_init <- rise
  run_init <- run
  
  # standardize the slope 
  rise <- rise / max(rise_init, run_init)
  run <- run / max(rise_init, run_init)
  
  # length that each projection traverses across each pixel
  l <- sqrt(rise^2 + run^2)
  
  #Return projection list
  #with d, length of projection within each pixel, indices (in order) beam went through, counts
  proj.list <- vector("list", 1)
  
  pl_idx <- 1
  
  
  
  # 1. Projections with slope (rise/run), starting from left wall
  for(rr in ROW){
    if(rr > 0){
      
      curr_row = rr
      curr_col = 1
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= floor(curr_col) && floor(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), floor(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), floor(curr_col)])
        curr_row = curr_row + rise
        curr_col = curr_col + run
        
      }
    }
    else{ # switch order of idx and theta if row number is negative
      curr_row = -rr
      curr_col = 1
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= floor(curr_col) && floor(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), floor(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), floor(curr_col)])
        curr_row = curr_row + rise
        curr_col = curr_col + run
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[!neg_ends(hit_thetas)] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[!neg_ends(hit_thetas)]
    
    if(length(idx) != 0 && sum(hit_thetas < 0) == 0){ 
      # if the projection didn't only go through -1's
      # and if there isn't a -1 in the middle of the vector
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  # 2. Projections with slope (-run/rise), starting from left wall
  for(rr in ROW){
    if(rr > 0){
      
      curr_row = rr
      curr_col = 1
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= ceiling(curr_col) && ceiling(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), ceiling(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), ceiling(curr_col)])    
        curr_row = curr_row - run
        curr_col = curr_col + rise
        
      }
    }
    else{ # switch order of idx and theta if row number is negative
      
      curr_row = -rr
      curr_col = 1
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= ceiling(curr_col) && ceiling(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), ceiling(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), ceiling(curr_col)])    
        curr_row = curr_row - run
        curr_col = curr_col + rise
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[!neg_ends(hit_thetas)] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[!neg_ends(hit_thetas)]
    
    if(length(idx) != 0 && sum(hit_thetas < 0) == 0){ 
      # if the projection didn't only go through -1's
      # and if there isn't a -1 in the middle of the vector
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  # 3. Projections with slope (rise/run), starting from bottom wall
  for(cc in COL[abs(COL) != 1]){ # skipping one to avoid redundancy
    if(cc > 0){
      
      curr_row = 1
      curr_col = cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= floor(curr_col) && floor(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), floor(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), floor(curr_col)])
        curr_row = curr_row + rise
        curr_col = curr_col + run
        
      }
    }
    else{ # switch order of idx and theta if col number is negative
      
      curr_row = 1
      curr_col = -cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= floor(curr_col) && floor(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), floor(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), floor(curr_col)])
        curr_row = curr_row + rise
        curr_col = curr_col + run
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[!neg_ends(hit_thetas)] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[!neg_ends(hit_thetas)]
    
    if(length(idx) != 0 && sum(hit_thetas < 0) == 0){ 
      # if the projection didn't only go through -1's
      # and if there isn't a -1 in the middle of the vector
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  # 4. Projections with slope (-run/rise), starting from top wall
  for(cc in COL[abs(COL) != 1]){ # skipping one to avoid redundancy
    if(cc > 0){
      
      curr_row = r
      curr_col = cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= ceiling(curr_col) && ceiling(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), ceiling(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), ceiling(curr_col)])
        curr_row = curr_row - run
        curr_col = curr_col + rise
        
      }
    }
    else{ # switch order of idx and theta if col number is negative
      
      curr_row = r
      curr_col = -cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= floor(curr_row) && floor(curr_row) <= r) && 
            (1 <= ceiling(curr_col) && ceiling(curr_col) <= c)){
        
        idx = c(idx, calc_index(r, floor(curr_row), ceiling(curr_col)))
        hit_thetas = c(hit_thetas, THETA[floor(curr_row), ceiling(curr_col)])
        curr_row = curr_row - run
        curr_col = curr_col + rise
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[!neg_ends(hit_thetas)] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[!neg_ends(hit_thetas)]
    
    if(length(idx) != 0 && sum(hit_thetas < 0) == 0){ 
      # if the projection didn't only go through -1's
      # and if there isn't a -1 in the middle of the vector
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  return(proj.list)
  
}





