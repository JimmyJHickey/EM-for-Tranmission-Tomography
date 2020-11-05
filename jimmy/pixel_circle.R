in_circle = function(radius){
  center = (2 * radius + 1 + 1) / 2
  nrow =  2 * radius + 1
  ncol =  2 * radius + 1
  
  # start everything at -1
  in_circ = matrix(-1, nrow = nrow, ncol = ncol)
  
  for(x in 1:nrow){
    for(y in 1:ncol){
      # if within the circle radius
      if( sqrt( (x - center)^2 + (y - center)^2 ) <= radius ){
        in_circ[x,y] = 0
      }
    }
  }
  
  # if in the circle, the value is 0, otherwise the value is -1
  return(in_circ)
}

