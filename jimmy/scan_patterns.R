rand_pattern = function(in_mat){
  out_mat = matrix(0, nrow = nrow(in_mat), ncol = ncol(in_mat))
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
     if(in_circ[x,y] == 1){
       out_mat[x,y] = abs(rnorm(1, 1, 2.5))
     } 
    }
  }
  
  return(out_mat)
}


line_pattern = function(in_mat, row = TRUE, line_number = 0){
  out_mat = matrix(0, nrow = nrow(in_mat), ncol = ncol(in_mat))

  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      if((in_circ[x,y] == 1) &&
         ((row && y == line_number) || (!row && x == line_number)) ){
        out_mat[x,y] = abs(rnorm(1, 1, 2.5))
      } 
    }
  }
  
  return(out_mat)
}

one_nonzero_pattern = function(in_mat, X, Y){
  out_mat = matrix(0, nrow = nrow(in_mat), ncol = ncol(in_mat))
  
  exit = FALSE
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      if((in_circ[x,y] == 1) &&
         (x == X && y == Y) ){
        out_mat[x,y] = abs(rnorm(1, 1, 2.5))
        exit = TRUE
        break
      } 
    }
    if(exit) break
  }
  
  return(out_mat)
}
