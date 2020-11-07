
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

#Set the names for the array
names.radius = c("rad3", "rad5", "rad10")
names.reps = c("reps1", "reps2", "reps3", "reps4", "reps5")
names.met = c("RMSE", "Spectral", "Iterations")
names.N = c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10")

names.list = list(names.radius, names.reps, names.met, names.N)

#Array to hold all of our results
#Rename to match your circle_theta name
checker_results <- array(NaN, dim =  c(3,5,3,10), dimnames = names.list)

#Output is a mulit-dimensional array
#Dim 1: Radius. Length=3 for radius=c(3,5,10)
#Dim 2: Reps. Length=5 for reps=1:5
#Dim 3: Metrics. Length=3 for metric = c('RMSE', 'Spectral', 'Iterations')
#Dim 4: Monte Carlo N. length=N=10
#Thus, we willl have an 3x5x3x10 array (=450 data points)

#Set the radius and rep numbers to loop over
radius.seq = c(3,5,10)
reps.seq = 1:5

##--------------------------------------------------------------
##--------------------------------------------------------------
for(radius in radius.seq){
  #Counter to know where to put the answers
  a = 1
  #Generating the circle_theta matrix. Different for each of us
  in_circ = in_circle(radius)
  circle_theta  = checker_pattern(in_circ)
  
  for(reps in reps.seq){
    b = 1
    for(N in 1:10){
      # generating observations
      set.seed(as.numeric(ceiling(proc.time()[3])))
      # vector of boundary rows/columns that aren't solely negative space
      # this code will not be right if there's more than one layer of -1's
      bounds <- 2:(nrow(circle_theta) - 1)
      
      #Generate data
      proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                               COL = c(bounds, -bounds),
                               reps = reps)
      
      #Check for no 0 values. If so, regenerate data
      while(y_zero(proj_list)){
        proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                                 COL = c(bounds, -bounds),
                                 reps = reps)        
      }
      
      # how many nonnegative numbers are in circle_theta?
      num_pixel <- sum(circle_theta >= 0)
      
      # copy the true circle_theta's negative space, but change the nonnegative
      # values randomly
      circle_theta_init <- circle_theta
      
      circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)
      
      #Run Algorithm
      em_res <- em_alg(proj_list, circle_theta_init, .0001) 
      

      
      abs_diff_mat <- abs(em_res$theta_est - circle_theta)
    
      sq_diff_mat <- abs_diff_mat^2
      
      #plot_matrix(abs_diff_mat)
      
      #plot_matrix(sq_diff_mat)
      
      # RMSE:
      checker_results[a,b,1,N] <- sqrt(sum(sq_diff_mat) / num_pixel) 
      #Spectral Norm
      checker_results[a,b,2,N] <-svd(em_res$theta_est - circle_theta)$d[1]
      #Number of iterations
      checker_results[a,b,3,N] <-em_res$ctr 
    }
    b = b+1  
  }
  a = a+1
}
  











