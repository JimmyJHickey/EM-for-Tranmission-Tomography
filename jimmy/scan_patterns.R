rand_pattern = function(in_mat){
  out_mat = in_mat
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
    # if in the circle
     if(in_mat[x,y] > 0){
       out_mat[x,y] = out_mat[x,y] + abs(rnorm(1, 1, 1))
     } 
    }
  }
  
  return(out_mat)
}


line_pattern = function(in_mat, row = TRUE, line_number = 0){
  out_mat = in_mat
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      # if in the circle and the correct row or column
      if((in_mat[x,y] > 0) &&
       ((row && y == line_number) || (!row && x == line_number)) ){
        out_mat[x,y] = out_mat[x,y] + abs(rnorm(1, 1, 2.5))
      } 
    }
  }
  
  return(out_mat)
}


one_nonzero_pattern = function(in_mat, X, Y){
  out_mat = in_mat
  
  # if in the circle 
  if(in_mat[X,Y] > 0){
    out_mat[X,Y] = out_mat[X,Y] + abs(rnorm(1, 1, 2.5))
  }
  
  return(out_mat)
}


checker_pattern = function(in_mat, row = TRUE, line_number = 0){
  out_mat = in_mat
  
  on = TRUE
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      # alternate on and off each pixel
      on = !on
      
      if((in_mat[x,y] > 0) ){
        # if on, generate a new data point
        if(on){
          out_mat[x,y] = out_mat[x,y] + abs(rnorm(1, 1, 1))
        }
      } 
    }
  }
  
  return(out_mat)
}


circle_pattern = function(in_mat, radius, center_x, center_y){
  out_mat = in_mat
  
  for(x in 1:nrow(out_mat)){
    for(y in 1:ncol(out_mat)){
      # if in the outer circle and the inscribed circle
      if((in_mat[x,y] >0) &&
        (sqrt( (x - center_x)^2 + (y - center_y)^2 ) <= radius )){
          # generate data inversely proportional to 
          # distance from center of inscribed circle
           out_mat[x,y] = out_mat[x,y] +  abs(rnorm(1, 
                              1/sqrt( 1+(x - center_x)^2 + (y - center_y)^2), 
                              1))
      } 
    }
  }
  
  return(out_mat)
}

