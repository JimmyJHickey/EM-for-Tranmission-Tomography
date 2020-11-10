load("em_results/boos/boos_results.RData")

results = boos_results

angles_vec = c(0,1,3)
rad_vec = c(10)

for(rad in 1:length(rad_vec)){
  cat(paste(rad_vec[rad], " & ", sep = ""))

  for (angle in 1:3){
    cat(paste(angles_vec[angle], " & ", sep=""))
        
    for(metric in 1:3){
      mean = round(mean(results[rad,angle,metric,]),3)
      stdev = round(var(sqrt(results[rad,angle,metric,])), 3)
      cat(paste(mean, " (", stdev,") & ", sep=""))
    }
    cat(" \\\\ \n")
  }
  cat("\\hline")
  
}


