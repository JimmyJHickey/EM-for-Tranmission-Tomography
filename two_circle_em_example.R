source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")



radius = 7

set.seed(1217)

in_circ = in_circle(radius)
plot_matrix(in_circ)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 2, 6, 6)
circle_theta2  = circle_pattern(circle_theta1, 2, 14, 14)
plot_matrix(circle_theta2)



circle_theta <- circle_theta2

reps <- 1

set.seed(1221)
# vector of boundary rows/columns that aren't solely negative space
# this code will not be right if there's more than one layer of -1's
bounds <- 2:(nrow(circle_theta) - 1)






source("Simulate Data.R")

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

save(em_res, file = "em_results/em_res.RData")

plot_matrix(em_res$theta_est)

abs_diff_mat <- abs(em_res$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

sqrt(sum(sq_diff_mat) / num_pixel) 
# 0 diag  0.08450741

# 1 diag  0.07840016

# 2 diag  0.06785882

# 4 diag  0.06785882

# 6 diag  0.06616761

# 8 diag  0.06921388

 
svd(em_r10_circle_reps5$theta_est - circle_theta)$d[1]
# 6.958846e+00







