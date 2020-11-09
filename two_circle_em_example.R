source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")

# The original function generates theta values that are too large. It's blocking all the photons in some projections.
# I modify the function to decrease the theta.

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
                                                 1)) # let's try shrinking the standard deviation
      } 
    }
  }
  
  return(out_mat)
}

radius = 7

set.seed(1217)

in_circ = in_circle(radius)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 2, 6, 6)
circle_theta2  = circle_pattern(circle_theta1, 2, 14, 14)
plot_matrix(circle_theta2)

circle_theta <- circle_theta2

reps <- 1





# first test, with 45 degree projections along with 90 degrees

set.seed(1221)

bounds <- 1:nrow(circle_theta)

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

length(proj_list) # 160

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res, file = "em_results/em_res.RData")

plot_matrix(em_res$theta_est)

abs_diff_mat <- abs(em_res$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) 
# 0.2767694

 
svd(em_res$theta_est - circle_theta)$d[1]
# 2.345076





# Trying things without 45 degree angles

set.seed(855)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds), rise_vec = c(), run_vec = c())

y_zero(proj_list)

length(proj_list) # 68

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_no_angle <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_no_angle, file = "em_results/em_res_no_angle.RData")

plot_matrix(em_res_no_angle$theta_est)

abs_diff_mat <- abs(em_res_no_angle$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) # the RMSE increases with no angles
# 0.4459833


svd(em_res_no_angle$theta_est - circle_theta)$d[1]
# 4.806766





# Trying things with projections with slopes of 1, 2 (along with perpendiculars), along with 90 degree angles

set.seed(910)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds), rise_vec = c(1, 2), run_vec = c(1, 1))

y_zero(proj_list)

length(proj_list) # 238

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_two_angles <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_two_angles, file = "em_results/em_res_two_angles.RData")

plot_matrix(em_res_two_angles$theta_est)

abs_diff_mat <- abs(em_res_two_angles$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) # the MSE decreases
# 0.1555816


svd(em_res_two_angles$theta_est - circle_theta)$d[1]
# 1.549698





# Trying things with projections with slopes of 1, 0.5 (along with perpendiculars), along with 90 degree angles

set.seed(930)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds), rise_vec = c(1, 2, 1), run_vec = c(1, 1, 2))

y_zero(proj_list)

length(proj_list) # 316

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_three_angles <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_three_angles, file = "em_results/em_res_three_angles.RData")

plot_matrix(em_res_three_angles$theta_est)

abs_diff_mat <- abs(em_res_three_angles$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) 
# 0.1183721


svd(em_res_three_angles$theta_est - circle_theta)$d[1]
# 1.164173




