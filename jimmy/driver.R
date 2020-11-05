radius = 20

in_circ = in_circle(radius)
plot_matrix(in_circ)

random_theta = rand_pattern(in_circ)
plot_matrix(random_theta)

line_theta = line_pattern(in_circ, FALSE, 5)
plot_matrix(line_theta)

one_theta = one_nonzero_pattern(in_circ, 10, 30)
plot_matrix(one_theta)

checker_theta  = checker_pattern(in_circ)
plot_matrix(checker_theta)

circle_theta  = circle_pattern(in_circ, 6, 25, 30)
plot_matrix(circle_theta)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 6, 20, 30)
circle_theta2  = circle_pattern(in_circ, 3, 30, 25)
plot_matrix(random_theta + circle_theta1 + circle_theta2)

