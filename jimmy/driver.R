radius = 10

in_circ = in_circle(radius)
plot_matrix(in_circ)

random_theta = rand_pattern(in_circ)
plot_matrix(random_theta)

line_theta = line_pattern(in_circ, FALSE, 5)
plot_matrix(line_theta)

one_theta = one_nonzero_pattern(in_circ, 4, 10)
plot_matrix(one_theta)

checker_theta  = checker_pattern(in_circ)
plot_matrix(checker_theta)

circle_theta  = circle_pattern(in_circ, 6, 25, 30)
plot_matrix(circle_theta)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 6, 4, 10)
circle_theta2  = circle_pattern(circle_theta1, 3, 30, 25)
plot_matrix(circle_theta2)



radius = 3

in_circ = in_circle(radius)
plot_matrix(in_circ)

# add two circles 
circle_theta1  = circle_pattern(in_circ, 2, 8, 8)
circle_theta2  = circle_pattern(circle_theta1, 2, 4, 4)
plot_matrix(circle_theta2)


radius = 3

in_circ = in_circle(radius)
plot_matrix(in_circ)

one_theta = one_nonzero_pattern(in_circ, 6, 4)
test = one_theta+one_theta+one_theta+one_theta+one_theta
plot_matrix(test)
