blank_plot <- function(text) {
  plot <- ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_text(aes(x = 0.5, y = 0.5, label = text), size = 6, hjust = 0.5, vjust = 0.5)+
    theme_bw()
  return(plot)
}