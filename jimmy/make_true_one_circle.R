set.seed(1979)

### Make faces

# radius 3
radius = 3

in_circ = in_circle(radius)
plot_matrix(in_circ)

one_circle = circle_pattern(in_circ, 1, 6, 7)
plot_matrix(one_circle)

save(one_circle, file = "true_theta/one_circle_rad3_truetheta.RData")
