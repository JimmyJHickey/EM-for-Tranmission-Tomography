set.seed(27605) # set the seed 
radius=10
in_circ = in_circle(radius)
circle_theta  = checker_pattern(in_circ)
plot_matrix(circle_theta)

b=1

rise.list = list(c(), c(1), c(1,2,1))
run.list = list(c(), c(1), c(1,1,2))  

set.seed(888)
# vector of boundary rows/columns that aren't solely negative space
# this code will not be right if there's more than one layer of -1's
bounds <- 1:nrow(circle_theta)

#Generate data
proj_list <- data_gen_df(circle_theta, d = 1000000000, ROW = c(bounds), 
                         COL = c(bounds),
                         rise_vec = rise.list[[b]], run_vec = run.list[[b]])
y_zero(proj_list)

# how many nonnegative numbers are in circle_theta?
num_pixel <- sum(circle_theta >= 0)

# copy the true circle_theta's negative space, but change the nonnegative
# values randomly
circle_theta_init <- circle_theta

circle_theta_init[which(circle_theta >= 0)] <- runif(num_pixel, 0, 0.1)

#Run Algorithm
em_res <- em_alg(proj_list, circle_theta_init, .0001) 
