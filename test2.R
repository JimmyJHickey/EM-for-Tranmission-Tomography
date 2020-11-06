
# read in all the functions

source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")

source("Simulate Data.R")

source("em_alg.R")



# generating the theta matrix

set.seed(627)

radius = 5

in_circ = in_circle(radius)
plot_matrix(in_circ)

circle_theta  = circle_pattern(in_circ, 2, 10, 10)
plot_matrix(circle_theta)



# generating observations

set.seed(631)

# data_gen_df returns numeric(0) if row/col is all -1. Automatic way to detect such row/cols?
# For this one, the nontrivial rows/cols are 2:14

# the EM algorithm breaks if any y are 0. I need to increase d until all ys are nonzero
# 1000000000 for the Poisson parameter works. Any larger, I get warnings.

proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(2:14, -c(2:14)), COL = c(2:14, -c(2:14)))



# running EM algorithm

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
(num_pixel <- sum(circle_theta >= 0))

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

em_r5_circle <- em_alg(proj_list, theta_init, .001) # less than 2 minutes

save(em_r5_circle, file = "em_results/em_r5_circle.RData")

plot_matrix(em_r5_circle$theta_est)

abs_diff_mat <- abs(em_r5_circle$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

mean(sq_diff_mat) # 0.2428365




### second test, reps = 2, otherwise same parameters as last test

# generating the theta matrix

set.seed(855)

radius = 5

in_circ = in_circle(radius)
# plot_matrix(in_circ)

circle_theta  = circle_pattern(in_circ, 2, 10, 10)
plot_matrix(circle_theta)

# generating observations

set.seed(900)

proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(2:14, -c(2:14)), COL = c(2:14, -c(2:14)),
                         reps = 2)

# running EM algorithm

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
(num_pixel <- sum(circle_theta >= 0))

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

system.time(em_r5_circle_reps2 <- em_alg(proj_list, theta_init, .001)) # took 1.03 minutes

save(em_r5_circle_reps2, file = "em_results/em_r5_circle_reps2.RData")

plot_matrix(em_r5_circle_reps2$theta_est)

abs_diff_mat <- abs(em_r5_circle_reps2$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

mean(sq_diff_mat) # 0.1727511

# adding one more repetition does decrease MSE





### third test, double the radius, same parameters as first test otherwise

# read in all the functions

source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")

source("Simulate Data.R")

source("em_alg.R")

# generating the theta matrix

set.seed(909)

radius = 10

in_circ = in_circle(radius)
plot_matrix(in_circ)

circle_theta  = circle_pattern(in_circ, 3, 10, 10)
plot_matrix(circle_theta)

# generating observations

set.seed(914)

# vector of boundary rows/columns that aren't solely negative space
bounds <- 2:(nrow(circle_theta) - 1)

proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds),
                         reps = 1)

# running EM algorithm

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
(num_pixel <- sum(circle_theta >= 0))

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

system.time(em_r10_circle <- em_alg(proj_list, theta_init, .0001)) # 16.39692 minutes

save(em_r10_circle, file = "em_results/em_r10_circle.RData")

plot_matrix(em_r10_circle$theta_est)

abs_diff_mat <- abs(em_r10_circle$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

mean(sq_diff_mat) # 0.1326812





### fourth test, 10x the repetitions, same parameters as last test otherwise

# read in all the functions

source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")

source("Simulate Data.R")

source("em_alg.R")

# generating the theta matrix

set.seed(909) # set the seed so that exactly same matrix is generated as the last test

radius = 10

in_circ = in_circle(radius)
plot_matrix(in_circ)

circle_theta  = circle_pattern(in_circ, 3, 10, 10)
plot_matrix(circle_theta)

# generating observations

set.seed(937)

# vector of boundary rows/columns that aren't solely negative space
bounds <- 2:(nrow(circle_theta) - 1)

proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds),
                         reps = 10)

# running EM algorithm

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
(num_pixel <- sum(circle_theta >= 0))

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

system.time(em_r10_circle_reps10 <- em_alg(proj_list, theta_init, .0001)) # 136.7684 minutes

save(em_r10_circle_reps10, file = "em_results/em_r10_circle_reps10.RData")

plot_matrix(em_r10_circle_reps10$theta_est)

abs_diff_mat <- abs(em_r10_circle_reps10$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

mean(sq_diff_mat) # 0.1909737




