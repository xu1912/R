library(plotly)
# Create a 3D scatter plot
plot_ly(data = d, x = ~PC1, y =~PC2, z = ~PC3, type = "scatter3d", mode = "markers", text=~Name)

t <- list(
  family = "sans serif",
  size = 14,
  color = toRGB("grey50"))

fig <- plot_ly(d, x = ~PC1, y =~PC2, z = ~PC3, text = ~Name)
fig <- fig %>% add_markers()
fig <- fig %>% add_text(textfont = t, textposition = "left")
fig <- fig %>%
    layout(
        xaxis = list(tickfont = list(size = 14)), 
        yaxis = list(tickfont = list(size = 14)),
	  zaxis = list(tickfont = list(size = 14)))

fig
