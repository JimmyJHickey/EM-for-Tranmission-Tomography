source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")

radius = 3

set.seed(238)

in_circ = in_circle(radius)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 1, 4, 8)
circle_theta2  = circle_pattern(circle_theta1, 1, 8, 8)
true_theta  = circle_pattern(circle_theta2, 1, 6, 4)
plot_matrix(true_theta)

# check the maximum thetas
summary(as.vector(true_theta))

save(true_theta, file = "true_theta/three_circles_rad3.RData")



# bounds <- 1:nrow(true_theta)
# 
# d <- 1e9
# 
# proj_list <- data_gen_df(true_theta, d = d, ROW = bounds,
#                          COL = bounds,
#                          rise_vec = c(1,2,1), run_vec = c(1,1,2))
# 
# 
# 
# # how many nonnegative numbers are in true_theta?
# num_pixel <- sum(true_theta >= 0)
# 
# # copy the true true_theta's negative space, but change the nonnegative
# # values randomly
# true_theta_init <- true_theta
# 
# true_theta_init[which(true_theta >= 0)] <- runif(num_pixel, 0, 0.1)
# 
# #Run Algorithm
# system.time(em_res <- em_alg(proj_list, true_theta_init, .0001))









