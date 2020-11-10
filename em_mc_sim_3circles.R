#PLEASE READ 
#This is our code to run our MC Experiments
#Each of us will run this code for a DIFFERENT circle_theta
#Checkered == Eric; One Circle == Alvin; Two Circles == Jimmy
#Thus, when we are all done, we will have one array each, three arrays total, one per circle_theta

##--------------------------------------------------------------
##--------------------------------------------------------------

#Read in Jimmy's Files
source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")

##--------------------------------------------------------------
##--------------------------------------------------------------

#Set the names for the array
names.radius = c("rad3", "rad5", "rad10")
names.angles = c("none", "deg45", "deg.all")
names.met = c("RMSE", "Spectral", "Iterations")
names.N = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10")

names.list = list(names.radius, names.angles, names.met, names.N)

# store the seeds
seed_vec <- rep(NA, length(names.radius) * length(names.angles) * length(names.N))
loop_idx <- 1

#Array to hold all of our results
#Rename to match your circle_theta name
three_circles_results <- array(NaN, dim =  c(3,3,3,10), dimnames = names.list)

#Output is a mulit-dimensional array
#Dim 1: Radius. Length=3 for radius=c(3,5,10)
#Dim 2: Angles Length=3 for angles: none, 45 deg and all three
#Dim 3: Metrics. Length=3 for metric = c('RMSE', 'Spectral', 'Iterations')
#Dim 4: Monte Carlo N. length=N=10
#Thus, we willl have an 3x5x3x10 array (=450 data points)

#Set the radius and rep numbers to loop over
radius.seq = c(3,5,10)
rise.list = list(c(), c(1), c(1,2,1))
run.list = list(c(), c(1), c(1,1,2))
a = 1
d = 1e9
##--------------------------------------------------------------
##--------------------------------------------------------------
for(radius in radius.seq){
  #Keeps track of progress
  print(paste("radius=",radius))
  #Generating the three_circles_theta matrix. Different for each of us
  load(paste("true_theta/three_circles_rad",radius,".RData",sep=""))
  three_circles_theta  = true_theta
  
  # iterate over angles
  for(b in 1:3){
    print(paste("angles_ind =",b))    
    for(N in 1:10){
      print(paste("N=",N))
      # generating observations
      seed <- as.numeric(ceiling(proc.time()[3]))
      seed_vec[loop_idx] <- seed
      set.seed(seed)
      # vector of boundary rows/columns that aren't solely negative space
      # this code will not be right if there's more than one layer of -1's
      bounds <- 1:nrow(three_circles_theta)
      
      #Generate data
      proj_list <- data_gen_df(three_circles_theta, d = d, ROW = bounds, 
                               COL = bounds,
                               rise_vec = rise.list[[b]], run_vec = run.list[[b]])
      
      #Check for no 0 values. If so, regenerate data
      while(y_zero(proj_list)){
        proj_list <- data_gen_df(three_circles_theta, d = d, ROW = bounds, 
                                 COL = bounds,
                                 rise_vec = rise.list[[b]], run_vec = run.list[[b]])        
      }
      
      # how many nonnegative numbers are in three_circles_theta?
      num_pixel <- sum(three_circles_theta >= 0)
      
      # copy the true three_circles_theta's negative space, but change the nonnegative
      # values randomly
      three_circles_theta_init <- three_circles_theta
      
      three_circles_theta_init[which(three_circles_theta >= 0)] <- runif(num_pixel, 0, 0.1)
      
      #Run Algorithm
      em_res <- em_alg(proj_list, three_circles_theta_init, .0001) 
      
      save(em_res, file=paste("em_results/three_circles/three_circles_radius",radius,"_numangles", b,"_N", N ,".RData",sep=""))
      
      abs_diff_mat <- abs(em_res$theta_est - three_circles_theta)
      
      sq_diff_mat <- abs_diff_mat^2
      
      #plot_matrix(abs_diff_mat)
      
      #plot_matrix(sq_diff_mat)
      
      # RMSE:
      three_circles_results[a,b,1,N] <- sqrt(sum(sq_diff_mat) / num_pixel) 
      #Spectral Norm
      three_circles_results[a,b,2,N] <-svd(em_res$theta_est - three_circles_theta)$d[1]
      #Number of iterations
      three_circles_results[a,b,3,N] <-em_res$ctr 
      
      loop_idx <- loop_idx + 1
      
    }
 
  }
  a = a+1
  
  save(three_circles_results, file = "three_circles_results.RData")
  print(three_circles_results)
  
}




save(three_circles_results, file = "three_circles_results.RData")



save(seed_vec, file = "seed_vec_three_circles_results.RData")









