plot_matrix = function(in_mat){
  require(ggplot2)

  x = seq(1, nrow(in_mat))
  y = seq(1, ncol(in_mat))
  
  df = expand.grid(X=x, Y=y)
  df$Z = c(in_mat)

  # doesn't seem to care about those colors at all, but it works
  ggplot(df, aes(X, Y, fill= Z)) + 
    geom_tile()+
    scale_color_manual(values = c("(-Inf,0]" = "black",
                                    "(0,Inf]" = "grey50"))
}
