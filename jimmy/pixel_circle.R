in_circle = function(radius){
  
  # radius of scan 
  # approximately 2% bigger than the radius of the patient cross section
  out_radius = as.integer(radius * 1.05 +1)
  
  center = (2 * out_radius + 1 + 1) %/% 2
  
  nrow =  2 * out_radius + 3
  ncol =  2 * out_radius + 3
  
  # start everything at -1
  in_circ = matrix(-1, nrow = nrow, ncol = ncol)
  
  # if outside of the scan -1
  # if inside the scan but outside the patient 0
  # if inside the patient random(Uniform(0, 0.25))
  for(x in 1:nrow){
    for(y in 1:ncol){
      # if within the patient radius
      if( sqrt( (x - center -1 )^2 + (y - center-1)^2 ) <= radius ){
        in_circ[x,y] = runif(1, 0, 0.25)
      }
      else if( sqrt( (x - center -1 )^2 + (y - center-1)^2 ) <= out_radius ){
        in_circ[x,y] = 0
      }
    }
  }
  

  return(in_circ)
}

