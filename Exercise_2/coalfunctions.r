library(scrm)
library(ape)

# SIM.TREE
# call sim.tree function to simulate gene trees
## This is to simulate a coalescent tree using the generation (discrete time) approximation algorithm, with 20 sampled lineages.
## In this example, the population size is 1000
# INPUT 
#   sample - number of lineages sampled. For instance, if 2 diploid individuals, sample=4.
#   current - number with current effective size (in number of diploids)
#   ancestral - number with ancestral effective size (in number of diploids)
#   time - number with time in generations ago of size change from ancestral to current
#   nrep - number of independent loci simulated
# OUTPUT
#   object of scrm with trees and segregating sites
sim.tree <- function(sample = 4,
                     current = 1000, 
                     ancestral=1000, 
                     time=0, nrep=1) {
  
  Nref <- current
  relA <- ancestral/Nref
  tchange <- time/(4*Nref)
  theta <- 4*Nref*1
  
  sum_stats <- scrm(paste(sample, nrep, '-t', theta ,'-T -eN', tchange, relA))
  sum_stats
}

# DRAW.TREE
## This plots (on the screen) the coalescent tree, with labels
# INPUT 
#   tree - object output from sim.tree
#   nrep - number of loci, i.e. independent simulations
#   events - value with scaled time of how long ago an event happened
# OUTPUT
#   draw trees with time measured in current Ne units, i.e. 1 unit is 4*Ne generations
draw.tree <- function(tree, nrep=1, events=-10) {
  # read the trees with ape package read.tree function
  trees <- read.tree(text = paste0(unlist(tree$trees)))
  
  if(nrep==1) {
    plot(trees, no.margin = FALSE)
    # add axis to tree
    axisPhylo()  
    
    # get the coalescent times
    coal_int <- coalescent.intervals(trees)
    tmrca <- sum(coal_int$interval.length)
    
    # add line of events
    abline(v=tmrca-events, col=4)
  } else if(nrep>1) {
    for(i in 1:length(trees)) {
      # plot the trees
      plot(trees[i], no.margin = FALSE)
      # add axis to tree
      axisPhylo()  
      
      # get the coalescent times
      coal_int <- coalescent.intervals(trees[[i]])
      tmrca <- sum(coal_int$interval.length)
      
      # add line of events
      abline(v=tmrca-events, col=4)
    }
  }
}

# COALTIMES
## Outputs coalescent times
# INPUT 
#   tree - object output from sim.tree
#   nrep - number of loci, i.e. independent simulations
# OUTPUT
#   matrix with each row corresponding to a coalescent interval
#   and each column to a simulation. The 1st column has the
#   number of lineages in each interval. NOTE: time scaled 
#   in units of 4*currentNe.
coaltimes <- function(tree, nrep=1) {
  # read the trees with ape package read.tree function
  trees <- read.tree(text = paste0(unlist(tree$trees)))
  
  if(nrep==1) {
    tmp <- coalescent.intervals(trees)
    coal_int <- matrix(NA, nrow=tmp$interval.count, ncol=2)
    coal_int[,1] <- tmp$lineages
    # get the coalescent times
    coal_int[,2] <- coalescent.intervals(trees)$interval.length
  } else if(nrep>1) {
    tmp <- coalescent.intervals(trees[[1]])
    coal_int <- matrix(NA, nrow=tmp$interval.count, ncol=length(trees)+1)
    coal_int[,1] <- tmp$lineages
    for(i in 1:length(trees)) {
      # get the coalescent times
      coal_int[,i+1] <- coalescent.intervals(trees[[i]])$interval.length
    }
  }
  coal_int
}


# SIM.TREE.MUT
## Simulates trees and adds mutations
# INPUT 
#   sample - number of lineages sampled. For instance, if 2 diploid individuals, sample=4.
#   current - number with current effective size (in number of diploids)
#   ancestral - number with ancestral effective size (in number of diploids)
#   time - number with time in generations ago of size change from ancestral to current
#   nrep - number of independent loci simulated
#   mu - number with mutation rate per site
#   L - integer number with number of sites 
# OUTPUT
#   matrix with each row corresponding to a coalescent interval
#   and each column to a simulation. The 1st column has the
#   number of lineages in each interval. NOTE: time scaled 
#   in units of 4*currentNe.
sim.tree.mut <- function(sample = 4,
                         current = 1000, 
                         ancestral=1000, 
                         time=0, nrep=1, 
                         mu=1e-8, L=1e4) {
  
  Nref <- current
  relA <- ancestral/Nref
  tchange <- time/(4*Nref)
  theta <- 4*Nref*L*mu
  sum_stats <- scrm(paste(sample, nrep, '-t', theta ,'-T -eN', tchange, relA))
  sum_stats
}


# PLOT.HAPLOTYPES
## This plots (on the screen) the coalescent tree, with labels
# INPUT 
#   tree - object output from sim.tree
#   nrep - number of loci, i.e. independent simulations
#   L    - integer number with number of sites
# OUTPUT
#   draw trees with time measured in current Ne units, i.e. 1 unit is 4*Ne generations
plot.haplotypes <- function(tree, nrep=1) {
  # read the trees with ape package read.tree function
  trees <- read.tree(text = paste0(unlist(tree$trees)))
  
  if(nrep==1) {
    par(mfrow=c(1,2))
    # plot gene tree  
    plot(trees, no.margin = FALSE)
    # add axis to tree
    axisPhylo()  
    
    # plot the SNPs
    # get the haplotypes
    hap <- tree$seg_sites[[1]]
    label <- as.numeric(trees$tip.label)
    # tmp <- matrix(0,ncol=L, nrow=nrow(hap))
    # tmp[,round(as.numeric(colnames(hap))*L)]  <- hap[as.numeric(trees$tip.label),]
    # image(x=1:ncol(tmp), y=1:nrow(tmp), t(tmp), ylab="lineages (inds)", xlab="sites")
    image(x=1:ncol(hap), y=1:nrow(hap), t(hap[label,]), ylab="lineages (inds)", xlab="SNPs", axes=FALSE)
    
  } else if(nrep>1) {
    for(i in 1:nrep) {
      par(mfrow=c(1,2))
      # plot the trees
      plot(trees[[i]], no.margin = FALSE)
      # add axis to tree
      axisPhylo()  
      # plot the SNPs
      # get the haplotypes
      hap <- tree$seg_sites[[i]]
      label <- as.numeric(trees[[i]]$tip.label)
      # tmp <- matrix(0,ncol=L, nrow=nrow(hap))
      # tmp[,round(as.numeric(colnames(hap))*L)]  <- hap[as.numeric(trees$tip.label),]
      # image(x=1:ncol(tmp), y=1:nrow(tmp), t(tmp), ylab="lineages (inds)", xlab="sites")
      image(x=1:ncol(hap), y=1:nrow(hap), t(hap[label,]), ylab="lineages (inds)", xlab="SNPs", axes=FALSE)
    }
  }
}

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
}
