# ABC_PLOT_REJECTION
# Function to plot the joint distribution of param and summary statistics, and show the points accepted after performing the rejection step. It also plots the histogram of the prior and posterior distributions.    
# INPUT:
# - simS: vector of size nsim with the summary statistic for each simulation
# - param: vector of size nsim with the parameter value used for each simulation
# - dist_sumstat: vector of size nsim with the distance between simulated and observed data, computed for each simulation.  
# - tolerance: value of the tolerance defined as the proportion of sims with smaller distance to the observed data (e.g. 0.10 indicates that we accept the parameters of the 10% closest simulations to the observed data).   
# OUTPUT:
# - plot with the joint distribution and histograms of prior and posterior distributions
abc_plot_rejection <- function(simS, param, dist_sumstat, tolerance) {
  # define plot parameters
  par(mfrow=c(1,1), mar=c(5,5,3,3))
  zones <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  # plot the joint distribution (x-axis: sumstat, y-axis: param)
  plot(simS, param, xlab="sumstat S", pch=1, ylab="Param")
  # get the tolerance as a quantile of distance distribution
  tol_dst <- quantile(dist_sumstat, tolerance)
  # get the index of simulations that are closer to observed data
  closest_sims <- which(dist_sumstat<tol_dst)
  # add points that are accepted coloured in red
  points(simS[closest_sims], param[closest_sims], col="red", pch=16)
  legend("topright", c("prior","posterior accepted"), pch=c(1,16), col=c(1,"red"))
  # add line corresponding to observed data
  abline(v=obs_S, col=4, lwd=2)    
  # Plot the histograms of the sumstat and parameters
  par(mar=c(0,5,3,3))
  breaks <- seq(min(simS),max(simS), length.out=100)
  x1hist <- hist(simS[closest_sims], plot=F, breaks=breaks)
  x2hist <- hist(simS, plot=F, breaks=breaks)
  ylims <- c(0, max(x2hist$density,x1hist$density))
  # Plot the distribution a priori of the statistic and the ones that are accepted
  barplot(x2hist$density, axes=T, ylim=ylims, space=0, horiz=F, col=1, border=1, add=F)
  barplot(x1hist$density, axes=T, ylim=ylims, space=0, horiz=F, col="red", border="red", add=T)
  par(mar=c(5,0,3,3)) 
  # plot the prior and the posterior
  breaks <- seq(min(param),max(param), length.out=100)
  x1hist <- hist(param[closest_sims], plot=F, breaks=breaks)
  x2hist <- hist(param, plot=F, breaks=breaks)
  ylims <- c(0, max(x1hist$density, x2hist$density))
  barplot(x2hist$density, axes=T, xlim=ylims, space=0, horiz=TRUE, col=1, border=1, add=F)
  barplot(x1hist$density, axes=T, xlim=ylims, space=0, horiz=TRUE, col="red", border="red", add=T)
  par(mfrow=c(1,1), mar=c(5,5,5,5))
}
