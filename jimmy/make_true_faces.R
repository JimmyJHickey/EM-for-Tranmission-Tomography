set.seed(1979)

### Make faces

# radius 3
radius = 3

in_circ = in_circle(radius)
plot_matrix(in_circ)

circ1 = circle_pattern(in_circ, 1, 4, 7)
circ2 = circle_pattern(circ1, 1, 8, 7)

face = line_pattern(circ2, TRUE, 4) 
plot_matrix(face)
max(face)

save(face, file = paste("true_theta/face_rad",radius,"_truetheta.RData", sep=""))


# radius 5
radius = 5

in_circ = in_circle(radius)
plot_matrix(in_circ)

circ1 = circle_pattern(in_circ, 2, 5, 10)
plot_matrix(circ1)
circ2 = circle_pattern(circ1, 2, 11, 10)
plot_matrix(circ2)

face = line_pattern(circ2, TRUE, 5) 
plot_matrix(face)
max(face)

save(face, file = paste("true_theta/face_rad",radius,"_truetheta.RData", sep=""))



# radius 10
radius = 10

in_circ = in_circle(radius)
plot_matrix(in_circ)

circ1 = circle_pattern(in_circ, 3, 8, 17)
plot_matrix(circ1)
circ2 = circle_pattern(circ1, 3, 18, 17)
plot_matrix(circ2)
face = line_pattern(circ2, TRUE, 9) 
plot_matrix(face)
max(face)

save(face, file = paste("true_theta/face_rad",radius,"_truetheta.RData", sep=""))

