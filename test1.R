
set.seed(708)

theta <- matrix(runif(9, 0, 1), nrow = 3)

dat <- data_gen_df(theta, d = 100, ROW = 1:3, COL = 1:3)

proj1 <- list(d = 100, idx = c(1, 4, 7), y = 6)

proj2 <- list(d = 100, idx = c(2, 5, 8), y = 28)

proj3 <- list(d = 100, idx = c(3, 6, 9), y = 21)

proj4 <- list(d = 100, idx = c(1, 2, 3), y = 19)

proj5 <- list(d = 100, idx = c(4, 5, 6), y = 24)

proj6 <- list(d = 100, idx = c(7, 8, 9), y = 10)

proj_list <- list(proj1, proj2, proj3, proj4, proj5, proj6)



theta_init <- matrix(runif(9, 0, 1), nrow = 3)



em_alg(proj_list, theta_init, .001)

