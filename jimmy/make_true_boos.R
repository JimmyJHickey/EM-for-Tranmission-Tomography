set.seed(1979)

### Make boos

# radius 10
radius = 10

boos = in_circle(radius)

# B
boos = one_nonzero_pattern(boos, 6, 10)
boos = one_nonzero_pattern(boos, 6, 11)
boos = one_nonzero_pattern(boos, 6, 12)
boos = one_nonzero_pattern(boos, 6, 13)
boos = one_nonzero_pattern(boos, 6, 14)
boos = one_nonzero_pattern(boos, 6, 15)
boos = one_nonzero_pattern(boos, 6, 16)


boos = one_nonzero_pattern(boos, 7, 10)
boos = one_nonzero_pattern(boos, 7, 13)
boos = one_nonzero_pattern(boos, 7, 16)

boos = one_nonzero_pattern(boos, 8, 11)
boos = one_nonzero_pattern(boos, 8, 12)
boos = one_nonzero_pattern(boos, 8, 14)
boos = one_nonzero_pattern(boos, 8, 15)


# O
boos = one_nonzero_pattern(boos, 10, 10)
boos = one_nonzero_pattern(boos, 10, 11)
boos = one_nonzero_pattern(boos, 10, 12)
boos = one_nonzero_pattern(boos, 10, 13)
boos = one_nonzero_pattern(boos, 10, 14)
boos = one_nonzero_pattern(boos, 10, 15)
boos = one_nonzero_pattern(boos, 10, 16)

boos = one_nonzero_pattern(boos, 11, 10)
boos = one_nonzero_pattern(boos, 11, 16)

boos = one_nonzero_pattern(boos, 12, 10)
boos = one_nonzero_pattern(boos, 12, 11)
boos = one_nonzero_pattern(boos, 12, 12)
boos = one_nonzero_pattern(boos, 12, 13)
boos = one_nonzero_pattern(boos, 12, 14)
boos = one_nonzero_pattern(boos, 12, 15)
boos = one_nonzero_pattern(boos, 12, 16)


# O
boos = one_nonzero_pattern(boos, 14, 10)
boos = one_nonzero_pattern(boos, 14, 11)
boos = one_nonzero_pattern(boos, 14, 12)
boos = one_nonzero_pattern(boos, 14, 13)
boos = one_nonzero_pattern(boos, 14, 14)
boos = one_nonzero_pattern(boos, 14, 15)
boos = one_nonzero_pattern(boos, 14, 16)

boos = one_nonzero_pattern(boos, 15, 10)
boos = one_nonzero_pattern(boos, 15, 16)

boos = one_nonzero_pattern(boos, 16, 10)
boos = one_nonzero_pattern(boos, 16, 11)
boos = one_nonzero_pattern(boos, 16, 12)
boos = one_nonzero_pattern(boos, 16, 13)
boos = one_nonzero_pattern(boos, 16, 14)
boos = one_nonzero_pattern(boos, 16, 15)
boos = one_nonzero_pattern(boos, 16, 16)

# S
boos = one_nonzero_pattern(boos, 18, 10)
boos = one_nonzero_pattern(boos, 18, 13)
boos = one_nonzero_pattern(boos, 18, 14)
boos = one_nonzero_pattern(boos, 18, 15)
boos = one_nonzero_pattern(boos, 18, 16)

boos = one_nonzero_pattern(boos, 19, 10)
boos = one_nonzero_pattern(boos, 19, 13)
boos = one_nonzero_pattern(boos, 19, 16)

boos = one_nonzero_pattern(boos, 20, 10)
boos = one_nonzero_pattern(boos, 20, 11)
boos = one_nonzero_pattern(boos, 20, 12)
boos = one_nonzero_pattern(boos, 20, 13)
boos = one_nonzero_pattern(boos, 20, 16)

plot_matrix(boos)
max(boos)

save(boos, file = paste("true_theta/boos_rad",radius,"_truetheta.RData", sep=""))


