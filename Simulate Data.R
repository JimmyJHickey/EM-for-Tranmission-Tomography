#ST793 Final Project
#EM algorithm for CT Data
#Eric Yanchenko (Alvin Sheng and Jimmy Hickey) [Copyright]
#November 4, 2020

calc_index = function(nrow, x, y){
  return(nrow * (y-1) + x)
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


#Function which generates the observed data for entire theta matrix and arbitrary projections
data_gen_df <- function(THETA, d, ROW, COL, reps=1){
  #THETA: input theta matrix used to compute probs.
  #d: Poisson mean for initial Poisson generation
  #ROW: rows which we project on (vector) Enter negative row number to get opposite direction
  #COL: columns which we project on (vector) Enter negative col number to get opposite direction
  #reps: number of times to run the function
  
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
    
    ###########
    # put diagonal projection function here
    # remember to increase the pl_idx accordingly
    ###########
    
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row+1
        curr_col = curr_col+1
        
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row+1
        curr_col = curr_col+1
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[hit_thetas >= 0] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[hit_thetas >= 0]
    
    if(length(idx) != 0){ # if the projection didn't only go through -1's
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])    
        curr_row = curr_row - 1
        curr_col = curr_col + 1
        
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])    
        curr_row = curr_row - 1
        curr_col = curr_col + 1
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[hit_thetas >= 0] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[hit_thetas >= 0]
    
    if(length(idx) != 0){ # if the projection didn't only go through -1's
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  # 3. Projections with slope (rise/run), starting from bottom wall
  for(cc in COL){
    if(cc > 0){
      
      curr_row = 1
      curr_col = cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row + 1
        curr_col = curr_col - 1
        
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row + 1
        curr_col = curr_col - 1
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[hit_thetas >= 0] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[hit_thetas >= 0]
    
    if(length(idx) != 0){ # if the projection didn't only go through -1's
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  # 4. Projections with slope (-run/rise), starting from top wall
  for(cc in COL){
    if(cc > 0){
      
      curr_row = 1
      curr_col = cc
      
      #Matrix indices that this beam goes through
      idx = c()
      
      # thetas intersected
      hit_thetas = c()
      
      # while still in bounds of the scan
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row + 1
        curr_col = curr_col + 1
        
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
      while((1 <= curr_row && curr_row <= r) && 
            (1 <= curr_col && curr_col <= c)){
        
        idx = c(idx, calc_index(r, curr_row, curr_col))
        hit_thetas = c(hit_thetas, THETA[curr_row, curr_col])
        curr_row = curr_row + 1
        curr_col = curr_col + 1
        
      }
      
      idx <- rev(idx)
      hit_thetas <- rev(hit_thetas)
      
    }
    
    idx = idx[hit_thetas >= 0] #Drop indices with negative thetas
    hit_thetas <- hit_thetas[hit_thetas >= 0]
    
    if(length(idx) != 0){ # if the projection didn't only go through -1's
      y <- data_gen(hit_thetas, d, l = l)
      add.list <- list(d = d, l = l, idx = idx, y = y)
      proj.list[[pl_idx]] <- add.list
      pl_idx <- pl_idx + 1
    }
    
  }
  
  return(proj.list)
  
}





