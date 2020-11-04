require("plot.matrix")

radius = 5

in_circ = in_circle(radius)
plot_matrix(in_circ)

random_theta = rand_pattern(in_circ)
plot_matrix(random_theta)

line_theta = line_pattern(in_circ, FALSE, 5)
plot_matrix(line_theta)

one_theta = one_nonzero_pattern(in_circ, 10, 3)
plot(one_theta)
