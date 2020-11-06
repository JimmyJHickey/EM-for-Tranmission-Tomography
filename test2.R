
# generating the theta matrix

source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")

set.seed(627)

radius = 5

in_circ = in_circle(radius)
plot_matrix(in_circ)

circle_theta  = circle_pattern(in_circ, 2, 10, 10)
plot_matrix(circle_theta)



# generating observations

source("Simulate Data.R")

set.seed(631)

# data_gen_df returns numeric(0) if row/col is all -1. Automatic way to detect such row/cols?
# For this one, the nontrivial rows/cols are 2:14

# the EM algorithm breaks if any y are 0. I need to increase d until all ys are nonzero
# 1000000000 for the Poisson parameter works. Any larger, I get warnings.

proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(2:14, -c(2:14)), COL = c(2:14, -c(2:14)))



# running EM algorithm

source("em_alg.R")

# there are no projections with zero observations
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true theta's negative space, but change the nonnegative
# values randomly
theta_init <- circle_theta

theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

em_alg(proj_list, theta_init, .001)
