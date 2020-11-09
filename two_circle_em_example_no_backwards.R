source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")

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
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds), 
                         COL = c(bounds),
                         reps = reps)

y_zero(proj_list)

length(proj_list) # 80

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_nob <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_nob, file = "em_results/em_res_nob.RData")

em_res_nob$ctr # 152

plot_matrix(em_res_nob$theta_est)

abs_diff_mat <- abs(em_res_nob$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) 
# 0.1908343


svd(em_res_nob$theta_est - circle_theta)$d[1]
# 1.758732




# Trying things without 45 degree angles

set.seed(851)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds), 
                         COL = c(bounds), rise_vec = c(), run_vec = c())

y_zero(proj_list)

length(proj_list) # 34

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_no_angle_nob <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_no_angle_nob, file = "em_results/em_res_no_angle_nob.RData")

em_res_no_angle_nob$ctr # 80

plot_matrix(em_res_no_angle_nob$theta_est)

abs_diff_mat <- abs(em_res_no_angle_nob$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) # the RMSE increases with no angles
# 0.4156665


svd(em_res_no_angle_nob$theta_est - circle_theta)$d[1]
# 4.560247





# Trying things with projections with slopes of 1, 2 (along with perpendiculars), along with 90 degree angles

set.seed(906)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds), 
                         COL = c(bounds), rise_vec = c(1, 2), run_vec = c(1, 1))

y_zero(proj_list)

length(proj_list) # 119

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_two_angles_nob <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_two_angles_nob, file = "em_results/em_res_two_angles_nob.RData")

em_res_two_angles_nob$ctr # 104

plot_matrix(em_res_two_angles_nob$theta_est)

abs_diff_mat <- abs(em_res_two_angles_nob$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) # the MSE decreases
# 0.1350274

svd(em_res_two_angles_nob$theta_est - circle_theta)$d[1]
# 1.225565





# Trying things with projections with slopes of 1, 0.5 (along with perpendiculars), along with 90 degree angles

set.seed(915)

bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds), 
                         COL = c(bounds), rise_vec = c(1, 2, 1), run_vec = c(1, 1, 2))

y_zero(proj_list)

length(proj_list) # 158

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res_three_angles_nob <- em_alg(proj_list, circle_theta_init, .0001) 

save(em_res_three_angles_nob, file = "em_results/em_res_three_angles_nob.RData")

em_res_three_angles_nob$ctr # 112

plot_matrix(em_res_three_angles_nob$theta_est)

abs_diff_mat <- abs(em_res_three_angles_nob$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# RMSE:

sqrt(sum(sq_diff_mat) / num_pixel) 
# 0.09653279

svd(em_res_three_angles_nob$theta_est - circle_theta)$d[1]
# 0.7616619


