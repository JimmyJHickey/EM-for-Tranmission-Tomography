
# DO EVERYTHING 10 TIMES 

### fifth test, 5x the repetitions, same parameters as last test otherwise. 

# read in all the functions

source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")

source("Simulate Data.R")

source("em_alg.R")

# generating the theta matrix

set.seed(909) # set the seed so that exactly same matrix is generated as the last test

# VARY THIS
radius = 10

in_circ = in_circle(radius)
plot_matrix(in_circ)

# CHANGE TO DIFFERENT THETA PATTERNS
circle_theta  = circle_pattern(in_circ, 3, 10, 10)
plot_matrix(circle_theta)

# generating observations

set.seed(749)

# vector of boundary rows/columns that aren't solely negative space
bounds <- 2:(nrow(circle_theta) - 1)

# VARY THIS
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds, -bounds), 
                         COL = c(bounds, -bounds),
                         reps = 5)

# running EM algorithm

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
(num_pixel <- sum(circle_theta >= 0))

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

system.time(em_r10_circle_reps5 <- em_alg(proj_list, theta_init, .0001)) # 92.35897 minutes

em_r10_circle_reps5$ctr # 324 iterations

save(em_r10_circle_reps5, file = "em_results/em_r10_circle_reps5.RData")

plot_matrix(em_r10_circle_reps5$theta_est)

abs_diff_mat <- abs(em_r10_circle_reps5$theta_est - circle_theta)

sq_diff_mat <- abs_diff_mat^2

plot_matrix(abs_diff_mat)

plot_matrix(sq_diff_mat)

# MSE:

sqrt(sum(sq_diff_mat)) / num_pixel # 0.167909
# 0.09065642


svd(em_r10_circle_reps5$theta_est - circle_theta)$d[1]
# 6.958846e+00


