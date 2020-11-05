
# test for 2x2 matrix

set.seed(740)

theta <- matrix(runif(4, 0, 1), nrow = 2)

dat <- data_gen_df(theta, d = 100000, ROW = 1:2, COL = 1:2)

proj1 <- list(d = 100000, idx = c(1, 3), y = 24)

proj2 <- list(d = 100000, idx = c(2, 4), y = 32)

proj3 <- list(d = 100000, idx = c(1, 2), y = 48)

proj4 <- list(d = 100000, idx = c(3, 4), y = 16)

proj_list <- list(proj1, proj2, proj3, proj4)



theta_init <- matrix(runif(4, 0, 1), nrow = 2)



em_alg(proj_list, theta_init, .00001)






# test for 3x3 matrix

set.seed(708)

theta <- matrix(runif(9, 0, 1), nrow = 3)

dat <- data_gen_df(theta, d = 100000, ROW = 1:3, COL = 1:3)

proj1 <- list(d = 100000, idx = c(1, 4, 7), y = 6)

proj2 <- list(d = 100000, idx = c(2, 5, 8), y = 28)

proj3 <- list(d = 100000, idx = c(3, 6, 9), y = 21)

proj4 <- list(d = 100000, idx = c(1, 2, 3), y = 19)

proj5 <- list(d = 100000, idx = c(4, 5, 6), y = 24)

proj6 <- list(d = 100000, idx = c(7, 8, 9), y = 10)

proj_list <- list(proj1, proj2, proj3, proj4, proj5, proj6)



theta_init <- matrix(runif(9, 0, 1), nrow = 3)



em_alg(proj_list, theta_init, .00001)






# test for 3x3 matrix

set.seed(708)

theta <- matrix(runif(9, 0, 1), nrow = 3)

dat <- data_gen_df(theta, d = 100000, ROW = 1:3, COL = 1:3)

proj1 <- list(d = 100000, idx = c(1, 4, 7), y = 6)

proj2 <- list(d = 100000, idx = c(2, 5, 8), y = 28)

proj3 <- list(d = 100000, idx = c(3, 6, 9), y = 21)

proj4 <- list(d = 100000, idx = c(1, 2, 3), y = 19)

proj5 <- list(d = 100000, idx = c(4, 5, 6), y = 24)

proj6 <- list(d = 100000, idx = c(7, 8, 9), y = 10)

dat2 <- data_gen_df(theta, d = 100000, ROW = 1:3, COL = 1:3)

proj7 <- list(d = 100000, idx = c(1, 4, 7), y = 10)

proj8 <- list(d = 100000, idx = c(2, 5, 8), y = 23)

proj9 <- list(d = 100000, idx = c(3, 6, 9), y = 32)

proj10 <- list(d = 100000, idx = c(1, 2, 3), y = 26)

proj11 <- list(d = 100000, idx = c(4, 5, 6), y = 15)

proj12 <- list(d = 100000, idx = c(7, 8, 9), y = 12)



proj_list <- list(proj1, proj2, proj3, proj4, proj5, proj6,
                  proj7, proj8, proj9, proj10, proj11, proj12)



theta_init <- matrix(runif(9, 0, 1), nrow = 3)



em_alg(proj_list, theta_init, .00001)


