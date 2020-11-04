rand_pattern = function(in_mat){
  out_mat = in_mat
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
     if(in_circ[x,y] == 1){
       out_mat[x,y] = abs(rnorm(1, 1, 2.5))
     } 
    }
  }
  
  return(out_mat)
}


line_pattern = function(in_mat, row = TRUE, num = 0){
  out_mat = in_mat
  
  line_number = if(num == 0) nrow(out_mat) / 2 else num
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      if((in_circ[x,y] == 1) &&
         ((row && y == num) || (!row && x == num)) ){
        out_mat[x,y] = abs(rnorm(1, 1, 2.5))
      } 
    }
  }
  
  return(out_mat)
}
