
uniroot(log, interval = c(0.0000001, 10), extendInt = "upX")$root  

uniroot(function(x) {-exp(x) + 5}, interval = c(0.0000001, 10), extendInt = "downX")$root  

uniroot(function(x) {exp(x) - 1e-16}, interval = c(0.0000001, 10), extendInt = "upX")$root  

uniroot(function(x) {1}, interval = c(0.0000001, 10), extendInt = "upX")$root  



# running code in IssuesWithEric.R

# hypothesis: observations too small? what if I made all projections traversing pixel 365 pretty small?

p <- length(proj_list)

proj_ys <- sapply(proj_list, function(proj) {proj$y})



(has_38 <- sapply(proj_list, function(proj) {38 %in% proj$idx}))

proj_ys[has_38]


q_fun_j(0.001, proj_list, circle_theta, 38) # 12 and 24



nonneg_idx <- which(circle_theta >= 0)

for (j in nonneg_idx[!(nonneg_idx %in% c(38, 59, 60, 61, 62, 64, 65, 66, 67, 82, 83, 93, 94, 106))]) {
  
  print(j)
  
  # what are the bounds of theta?
  # if I need more efficiency, calculate nij and mij outside of function
  uniroot(q_fun_j, interval = c(0.0000001, 10), proj_list, circle_theta, j, 
                          extendInt = "downX")$root  
  
}


# projection 24
nij(proj_list[[24]], circle_theta, 38)
mij(proj_list[[24]], circle_theta, 38)

# projection 12
nij(proj_list[[12]], circle_theta, 38)
mij(proj_list[[12]], circle_theta, 38)




# pixel 65

(has_65 <- sapply(proj_list, function(proj) {65 %in% proj$idx}))

c(1:p)[has_65] # 14 and 25

proj_ys[has_65]



# projection 14
nij(proj_list[[14]], circle_theta, 65)
mij(proj_list[[14]], circle_theta, 65)

# projection 25
nij(proj_list[[25]], circle_theta, 65)
mij(proj_list[[25]], circle_theta, 65)




###########

# pixel that works

(has_110 <- sapply(proj_list, function(proj) {110 %in% proj$idx}))

c(1:p)[has_110] # 9 27

proj_ys[has_110]

# projection 1
nij(proj_list[[9]], circle_theta, 110)
mij(proj_list[[9]], circle_theta, 110)

# projection 35
nij(proj_list[[27]], circle_theta, 110)
mij(proj_list[[27]], circle_theta, 110)

#
uniroot(q_fun_j, interval = c(0.0000001, 10), proj_list, circle_theta, 110, 
        extendInt = "downX")$root  





