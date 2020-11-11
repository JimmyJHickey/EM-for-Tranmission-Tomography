
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

set.seed(27605) # set the seed 

NN = 100

#Set the names for the array
names.radius = c("rad3")
names.angles = c("none", "deg45", "deg.all")
names.met = c("RMSE", "Spectral", "Iterations")
names.N = seq(1:NN)

names.list = list(names.radius, names.angles, names.met, names.N)

#Array to hold all of our results
#Rename to match your circle_theta name
one_circle_results <- array(NaN, dim =  c(1,3,3,NN), dimnames = names.list)

#Output is a mulit-dimensional array
#Dim 1: Radius. Length=3 for radius=c(3,5,10)
#Dim 2: Angles Length=3 for angles: none, 45 deg and all three
#Dim 3: Metrics. Length=3 for metric = c('RMSE', 'Spectral', 'Iterations')
#Dim 4: Monte Carlo N. length=N=10
#Thus, we willl have an 3x5x3x10 array (=450 data points)

#Set the radius and rep numbers to loop over
radius.seq = c(3)
rise.list = list(c(), c(1), c(1,2,1))
run.list = list(c(), c(1), c(1,1,2))
a = 1
d = 1000000000
##--------------------------------------------------------------
##--------------------------------------------------------------
for(radius in radius.seq){
  #Keeps track of progress
  print(paste("radius=",radius))
  #Generating the one_circle_theta matrix. Different for each of us
  load(paste("true_theta/one_circle_rad",radius,"_truetheta.RData",sep=""))
  one_circle_theta  = one_circle
  
  # iterate over angles
  for(b in 1:3){
    print(paste("angles_ind =",b))    
    for(N in 1:NN){
      print(paste("N=",N))
      # generating observations
      set.seed(as.numeric(ceiling(proc.time()[3])))
      # vector of boundary rows/columns that aren't solely negative space
      # this code will not be right if there's more than one layer of -1's
      bounds <- 1:nrow(one_circle_theta)
      
      #Generate data
      proj_list <- data_gen_df(one_circle_theta, d = d, ROW = bounds, 
                               COL = bounds,
                               rise_vec = rise.list[[b]], run_vec = run.list[[b]])
      
      #Check for no 0 values. If so, regenerate data
      while(y_zero(proj_list)){
        proj_list <- data_gen_df(one_circle_theta, d = d, ROW = bounds, 
                                 COL = bounds,
                                 rise_vec = rise.list[[b]], run_vec = run.list[[b]])    
      }
      
      # how many nonnegative numbers are in one_circle_theta?
      num_pixel <- sum(one_circle_theta >= 0)
      
      # copy the true one_circle_theta's negative space, but change the nonnegative
      # values randomly
      one_circle_theta_init <- one_circle_theta
      
      one_circle_theta_init[which(one_circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)
      
      #Run Algorithm
      em_res <- em_alg(proj_list, one_circle_theta_init, .0001) 
      
      save(em_res, file=paste("em_results/one_circle/one_circle_radius",radius,"_numangles", b,"_N", N ,".RData",sep=""))
      
      abs_diff_mat <- abs(em_res$theta_est - one_circle_theta)
      
      sq_diff_mat <- abs_diff_mat^2
      
      #plot_matrix(abs_diff_mat)
      
      #plot_matrix(sq_diff_mat)
      
      # RMSE:
      one_circle_results[a,b,1,N] <- sqrt(sum(sq_diff_mat) / num_pixel) 
      #Spectral Norm
      one_circle_results[a,b,2,N] <-svd(em_res$theta_est - one_circle_theta)$d[1]
      #Number of iterations
      one_circle_results[a,b,3,N] <-em_res$ctr 
    }
    b = b+1  
  }
  a = a+1
}

save(one_circle_results, file = "em_results/one_circle/one_circle_results.RData")













