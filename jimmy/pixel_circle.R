library("plot.matrix")

radius = 5
center = (2 * radius + 1 + 1) / 2
nrow =  2 * radius + 1
ncol =  2 * radius + 1

in_circ = matrix(0, nrow = nrow, ncol = ncol)

for(i in 1:nrow){
  for(j in 1:ncol){
    if( sqrt( (i - center)^2 + (j - center)^2 ) <= radius ){
      in_circ[i,j] = 1
    }
  }
}

plot(in_circ)

