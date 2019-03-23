# @description Plot the dFC time dependence  modeled in dynamic functional connectivity modelling
# @title 95\% Coverage time dependence plots for dFCM
# @param data.p A matrix with pairwise correlation estimates at each time points and corresponding confidence interval generated from lme modeling.
# @param roi1 index of 1st region of interest, eg. 278
# @param roi2 index 2nd region of interest, eg. 274
# @param n_timepoints number of total time.points per scan
# @param lag shift unit between two pairwise time series, eg. "-1", "-2", "+1", "0" etc.
# @family plot
# 
#



plot_timeDependence_dFCM <- function(data.p,
                                     roi1,
                                     roi2,
                                     n_timepoints) {

  #set_plotter_package("curve") # Load required packages for plotting
  
  requireNamespace("ggplot2")


  ind1 <- which(roi_names[, 1] == as.integer(roi1))
  ind2 <- which(roi_names[, 1] == as.integer(roi2))

  title_plot <- paste(roi_names[ind1, 6], "vs", roi_names[ind2, 6])


  # Conigure x axis (timepoints)
  x <- seq(0.5, n_timepoints, length.out = ncol(data.p))


  plotly_data_p <- data.frame('x' = x,
                              'y_diff' = as.numeric(data.p[2, ]),
                              'y_beer' = as.numeric(data.p[1, ]),
                              'y_gat' = as.numeric(data.p[3, ]),
                              'y_diff_lowb' = as.numeric(data.p[4, ]),
                              'y_diff_upperb' = as.numeric(data.p[5, ]),
                              'y_beer_lowb' = as.numeric(data.p[6, ]),
                              'y_beer_upperb' = as.numeric(data.p[7, ]),
                              'y_gat_lowb' = as.numeric(data.p[8, ]),
                              'y_gat_upperb' = as.numeric(data.p[9, ]))
  #print(plotly_data_p)
  # Initialize an empty figure
  p <- ggplot2::ggplot(plotly_data_p) + ggplot2::geom_hline(yintercept = 0, size = I(0.8))  +  # Horizontal line at x = 0
    ggplot2::geom_line(ggplot2::aes(x = plotly_data_p$x, y = plotly_data_p$y_gat, colour = "Cond 2", size = I(0.5)), show.legend = T) +   # Line for condition 2
    ggplot2::geom_line(ggplot2::aes(x = plotly_data_p$x, y = plotly_data_p$y_beer, colour = "Cond 1", size = I(0.5)), show.legend = T) +  # Line for condition 1
    ggplot2::geom_line(ggplot2::aes(x = plotly_data_p$x, y = plotly_data_p$y_diff, colour = "Cond Diff.", size = I(0.5)), show.legend = T)+  # Line for condition2 - condition 1
    ggplot2::geom_ribbon(ggplot2::aes(ymin = plotly_data_p$y_beer_lowb, ymax = plotly_data_p$y_beer_upperb, x = x, fill = "Cond 1"), alpha = 0.2) +  # Shaded area
    ggplot2::geom_ribbon(ggplot2::aes(ymin = plotly_data_p$y_gat_lowb, ymax = plotly_data_p$y_gat_upperb, x = x, fill = "Cond 2"), alpha = 0.2) +   # Shaded area
    ggplot2::geom_ribbon(ggplot2::aes(ymin = plotly_data_p$y_diff_lowb, ymax = plotly_data_p$y_diff_upperb, x = x, fill = "Cond Diff."), alpha = 0.2) +   # Shaded area
    ggplot2::scale_x_continuous(limits = c(0.5,n_timepoints), expand = c(0,0))+ # Limit of x and y
    ggplot2::ggtitle(label = title_plot) + ggplot2::xlab("Time Points") + ggplot2::ylab("Corr. Estimate")+ # Information of the plot
    ggplot2::theme(plot.title = ggplot2::element_text(size = 20,
                                    lineheight = 1.2,
                                    face = "bold"),
          axis.title = ggplot2::element_text(size = 15,
                                   lineheight = 1.2),
          legend.text = ggplot2::element_text(size = 10,
                                     lineheight = 1.2),
          legend.title = ggplot2::element_text(size = 15,
                                     lineheight = 1.2))+
    ggplot2::scale_colour_manual(values = c("Cond 1" = "blue", "Cond 2" = "green", "Cond Diff." = "red"), aesthetics = c("colour", "fill"))


return(p)

}

# @importFrom("grDevices", "dev.off")
# @importFrom("graphics", "abline", "axis", "box", "text")
# @description Plot the pairwise non-zero coverage matrix modeled in dynamic functional connectivity modelling
# @title Non-zero coverage of Pairwise comparison.
# @param mat non-zero coverage matrix to plot.
# @param file filepath to generate plot.
# @param label.x label of x-axis
# @param label.y label of y-axis
# @param title title of the figure
#@param type type of plots, c("combined", "separate")
# @param height height of the figure(in inch), default is 17 inch
# @param width width of the figure(in inch), default is 17 inch
#@param res resolution of the figure
# @param mar margin of the figure.
# @importFrom grDevices dev.off 
# @importFrom grDevices png 
# @importFrom graphics axis abline
# @importFrom graphics box
# @importFrom graphics image
# @importFrom graphics par
# @importFrom graphics text
# @importFrom gplots colorpanel
# @family plot

plot_non_zero_coverage <- function(mat,
                                   save_fig,
                                   file,
                                   label.x,
                                   label.y,
                                   title,
                                   type = "combined",
                                   height = 17,
                                   width = 17,
                                   res = 700,
                                   mar = c(14, 14, 3.5, 6.5)){
  
  lp.roi <- nrow(mat)
  tab1 <- rep(1, lp.roi)
  tab2 <- cumsum(tab1)
  tab2 <- c(0, (tab2))
  tab2.tick <- seq(0, lp.roi, 0.5)
  
  if(isTRUE(save_fig)){
    grDevices::png(file,
                   width = width,
                   height = height,
                   units = "in",
                   res = res)
  }
  
  # Set up the margin
  graphics::par(mar = mar)
  
  if(type == "combined"){
    graphics::image(c(0:lp.roi),
          c(0:lp.roi),
          ylim = c(lp.roi, 0),
          as.matrix(mat),
          zlim = c(-0.01, 1),
          main = title,
          col = c("blue", "grey", gplots::colorpanel(100, "white", "salmon3", "gold")),
          axes = F,
          xlab = "",
          ylab = "")
    
  }else{
    graphics::image(c(0:lp.roi), c(0:lp.roi),
          ylim = c(lp.roi, 0),
          as.matrix(mat),
          zlim = c(-0.01, 1),
          main = title,
          col = c("blue", "grey", gplots::colorpanel(100, "white", "salmon3", "gold")),
          axes = F,
          xlab = "",
          ylab = "")
  }

  graphics::box()
  
  graphics::axis(2, at = tab2.tick, labels = label.y, cex.axis = 0.05*width, las = 2)
  graphics::axis(2, at = tab2.tick, labels = label.y, cex.axis = 0.05*width, col = "white", las = 2)
  graphics::axis(2, at = tab2, labels = rep(" ", length(tab2)), cex.axis = 0.05*width, las = 2)
  graphics::axis(1, at = tab2.tick, labels = label.x, cex.axis = 0.05*width, las = 2)
  graphics::axis(1, at = tab2.tick, labels = label.x, cex.axis = 0.05*width, col = "white", las = 2)
  graphics::axis(1, at = tab2, labels = rep(" ", length(tab2)), cex.axis = 0.05*width, las = 2)
  
  graphics::abline(v = (tab2), col = "black", lty = "dotted")
  graphics::abline(h = (tab2), col = "black", lty = "dotted")
  
  if(type=="combined"){
    fields::image.plot(as.matrix(mat),
               zlim = c(-0.01, 1),
               main = title,
               col = c("blue", "grey", gplots::colorpanel(100, "white", "salmon3", "gold")),
               legend.only = T)    
  }else{
    fields::image.plot(as.matrix(mat), zlim = c(-0.01, 1), main = "Non-zero coverage for condtion1 and condition2",
               col = c("blue", "grey", gplots::colorpanel(100, "white", "salmon3", "gold")), legend.only = T)
    
  }

  
  
  diag(mat) <- 0
  
  for (x in seq(0.5, lp.roi, 1)) for (y in seq(0.5, lp.roi, 1)) if (x != y) {
    graphics::text(x, y, t(as.matrix(round(mat, 2)))[(y + 0.5), (x + 0.5)], cex = 0.1*width)
  }
  
  
  
  graphics::box()
  if(save_fig) grDevices::dev.off()
}
