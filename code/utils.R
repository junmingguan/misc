# credit to Annie
plot_heatmap <- function(L, title = "", colors_range = c("blue", "white", "red"), brks = NULL){
  ### define the color map
  cols <- colorRampPalette(colors_range)(49)
  if (is.null(brks) == TRUE){
    lim <- max(abs(L), na.rm = TRUE)
    brks <- seq(-lim, lim, length=50)
  }

  plt <- pheatmap(L, show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE, color = cols, breaks = brks, main = title)
  return(plt)
}


harminize_L_sign = function(Lhat, obj) {
  ord <- order(obj, decreasing = FALSE)
  Lhat_ord <- Lhat[ord, ]

  ref <- Lhat_ord[1, ]

  sgn <- sign(drop(Lhat_ord %*% ref))
  sgn[sgn == 0] <- 1

  Lhat_ord * sgn
}
