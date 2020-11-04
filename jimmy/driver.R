require("plot.matrix")

radius = 5

in_circ = in_circle(radius)
plot_matrix(in_circ)

random_theta = rand_pattern(in_circ)
plot_matrix(random_theta)

line_theta = line_pattern(in_circ, FALSE, 3)
plot_matrix(line_theta)
