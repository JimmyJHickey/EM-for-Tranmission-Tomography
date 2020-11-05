
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

proj_list <- data_gen_df(circle_theta, d = 100000, ROW = c(2:14, -c(2:14)), COL = c(2:14, -c(2:14)))




