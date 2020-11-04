plot_matrix = function(in_mat){
  require(ggplot2)

  x = paste0("X", seq(1, nrow(in_mat)))
  y = paste0("Y", seq(1, ncol(in_mat)))
  
  df = expand.grid(X=x, Y=y)
  df$Z = c(in_mat)

  
  ggplot(df, aes(X, Y, fill= Z)) + 
    geom_tile()
}
