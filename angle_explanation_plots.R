source("jimmy/pixel_circle.R")
source("jimmy/plotting.R")
source("jimmy/scan_patterns.R")
source("Simulate Data.R")
source("em_alg.R")

library(ggplot2)



# adjusted function to emphasize inside of patient 

in_circle = function(radius){
  
  # radius of scan 
  # approximately 5% bigger than the radius of the patient cross section
  out_radius = as.integer(radius * 1.05 +1)
  
  center = (2 * out_radius + 1 + 1) %/% 2
  
  nrow =  2 * out_radius + 3
  ncol =  2 * out_radius + 3
  
  # start everything at -1
  in_circ = matrix(-1, nrow = nrow, ncol = ncol)
  
  # if outside of the scan -1
  # if inside the scan but outside the patient 0
  # if inside the patient random(Uniform(0, 0.25))
  for(x in 1:nrow){
    for(y in 1:ncol){
      # if within the patient radius
      if( sqrt( (x - center -1 )^2 + (y - center-1)^2 ) <= radius ){
        in_circ[x,y] = 4
      }
      else if( sqrt( (x - center -1 )^2 + (y - center-1)^2 ) <= out_radius ){
        in_circ[x,y] = 0
      }
    }
  }
  
  
  return(in_circ)
}



radius = 3

set.seed(1217)

in_circ = in_circle(radius)

plot_matrix(in_circ)



# modified to add lines

x = seq(1, nrow(in_circ))
y = seq(1, ncol(in_circ))

df = expand.grid(X=x, Y=y)
df$Z = c(in_circ)

# doesn't seem to care about those colors at all, but it works

ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_continuous(type = "viridis",
                        limits = c(0,4), 
                        breaks = c(0, 1, 2, 3, 4),
                        guide_colourbar(nbin = 100),
                        name = "theta") + 
  geom_hline(yintercept = 2:10, 
             color = "red", size = 0.5) + 
  geom_vline(xintercept = 2:10, 
             color = "red", size = 0.5)



ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_continuous(type = "viridis",
                        limits = c(0,4), 
                        breaks = c(0, 1, 2, 3, 4),
                        guide_colourbar(nbin = 100),
                        name = "theta") + 
  geom_hline(yintercept = 2:10, 
             color = "red", size = 0.5) + 
  geom_vline(xintercept = 2:10, 
             color = "red", size = 0.5) + 
  geom_abline(intercept = -4:4, slope = 1, color = "blue", size = 0.5) + 
  geom_abline(intercept = 8:16, slope = -1, color = "blue", size = 0.5)



ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_continuous(type = "viridis",
                        limits = c(0,4), 
                        breaks = c(0, 1, 2, 3, 4),
                        guide_colourbar(nbin = 100),
                        name = "theta") + 
  geom_hline(yintercept = 2:10, 
             color = "red", size = 0.5) + 
  geom_vline(xintercept = 2:10, 
             color = "red", size = 0.5) + 
  geom_abline(intercept = -5:5, slope = 1, color = "blue", size = 0.5) + 
  geom_abline(intercept = 8:16, slope = -1, color = "blue", size = 0.5) +
  geom_abline(intercept = -13:1, slope = 2, color = "green3", size = 0.5) +
  geom_abline(intercept = 6:12, slope = -1/2, color = "green3", size = 0.5) + 
  geom_abline(intercept = -13:1, slope = 1/2, color = "green3", size = 0.5) +
  geom_abline(intercept = 6:12, slope = -2, color = "green3", size = 0.5) + 
  geom_abline(intercept = -13:1, slope = 2, color = "green3", size = 0.5) +
  geom_abline(intercept = 6:12, slope = -1/2, color = "green3", size = 0.5) +
  geom_abline(intercept = 0:6, slope = 1/2, color = "green3", size = 0.5) +
  geom_abline(intercept = 11:25, slope = -2, color = "green3", size = 0.5)



ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile()+
  scale_fill_continuous(type = "viridis",
                        limits = c(0,4), 
                        breaks = c(0, 1, 2, 3, 4),
                        guide_colourbar(nbin = 100),
                        name = "theta") + 
  # geom_hline(yintercept = 2:10, 
  #            color = "red", size = 0.5) + 
  # geom_vline(xintercept = 2:10, 
  #            color = "red", size = 0.5) + 
  # geom_abline(intercept = -5:5, slope = 1, color = "blue", size = 0.5) + 
  # geom_abline(intercept = 8:16, slope = -1, color = "blue", size = 0.5) +
  geom_abline(intercept = -13:1, slope = 2, color = "green3", size = 0.5) +
  geom_abline(intercept = 6:12, slope = -1/2, color = "green3", size = 0.5) +
  geom_abline(intercept = 0:6, slope = 1/2, color = "green3", size = 0.5) +
  geom_abline(intercept = 11:25, slope = -2, color = "green3", size = 0.5)





abline(h = 1:11, col = "red", lwd = 3)
abline(v = 1:11, col = "red", lwd = 3)






# 4 x 4 example matrix

