plot_matrix = function(in_mat){
  require(ggplot2)

  x = seq(1, nrow(in_mat))
  y = seq(1, ncol(in_mat))
  
  df = expand.grid(X=x, Y=y)
  df$Z = c(in_mat)

  # doesn't seem to care about those colors at all, but it works
  ggplot(df, aes(X, Y, fill= Z)) + 
    geom_tile()+
    scale_fill_continuous(type = "viridis",
                          limits = c(0,4), 
                          breaks = c(0, 1, 2, 3, 4),
                          guide_colourbar(nbin = 100),
                          name = "theta")
  }
