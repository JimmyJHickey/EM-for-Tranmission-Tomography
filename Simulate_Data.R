#ST793 Final Project
#EM algorithm for CT Data
#Eric Yanchenko (Alvin Sheng and Jimmy Hickey) [Copyright]
#November 4, 2020

#Function which generates the observed data for a single theta vector
data_gen <- function(theta, d){
  #Theta is the given theta vector used to compute attenutation probs.
  #d is mean of initial Poisson distribution
  
  p = length(theta)
  
  X1 = rpois(1, lambda = d) #Initial number of counts from Poisson dist.
  
  X = c(X1, rep(0, p)) #Vector to hold the number of counts after each interaction
  
  for(j in 1:p){
    X[j+1] <- rbinom(1, size = X[j], prob = exp(-theta[j]))
    #Generate the number of photons from Binom(X_i, exp(theta_i))
  }
  
  Y = X[p+1] #Actually observed data (number of counts)
  
  return(Y)
}

#Function which generates the observed data for entire theta matrix and arbitrary projections
data_gen_df <- function(THETA, d, ROW, COL){
  #THETA: inputed theta matrix used to compute probs.
  #d: Poisson mean for initial Poisson generation
  #ROW: rows which we project on (vector)
  #COL: columns which we project on (vector)
  
  r = dim(THETA)[1] #Number of rows in THETA
  c = dim(THETA)[2] #Number of cols in THETA
  
  row.dat <- numeric(length(ROW)) #Holds the row data
  col.dat <- numeric(length(COL)) #Holds the col data
  
  #Check that ROW and COL projection indices are valid
  if(max(ROW) > r || max(COL) > c){
    return(simpleError("Row and/or column indices out of range."))
  }
  
  #Run over ROWS first
  i = 1
  for(rr in ROW){
    row.dat[i] <- data_gen(THETA[rr,], d)
    i = i+1
  }
  
  #Run over COLS second
  j = 1
  for(cc in COL){
    col.dat[j] <- data_gen(THETA[,cc], d)
    j = j+1
  }
  
  #Return the data for rows and columns as well as the indices and Poisson mean parameter
  return(list(row.dat = row.dat, row.ind =  ROW, col.dat = col.dat, col.ind = COL, lambda = d))
}







