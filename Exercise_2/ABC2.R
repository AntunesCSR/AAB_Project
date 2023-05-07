# AAB 2023
# Individual Assignment: Introduction to Approximate Bayesian Computation (ABC)
# Author: Cátia Antunes
# Date: 25-04-2023

##########################################################################################################################

# Install packages
install.packages("locfit")
install.packages("scrm")
install.packages("ape")
install.packages("scrm")

# Load the packages
library("locfit")
library("scrm")
library("ape")
library(scrm)
source("coalfunctions.r")
source("coalfunctions.r")
source("abc_plot_rejection.r")

##########################################################################################################################

setwd("./Exercise_2") # set working directory

### 1. Select data to be used in the analysis.

# Load data
pan_data <- read.table("SegregatingSitesChimp.txt", header = TRUE)

pop <- pan_data$SpeciesPop # load population
pop
ind_num <- pan_data$NumIndividuals[1] # load number of individuals
ind_num
seg_sites <- pan_data$SegregatingSites # load segregating sites
seg_sites
num_sites <- pan_data$NumSites[1] # load number of sites
num_sites
mt_rate <- pan_data$MtRatePerSite[1] # load mutation rate per site
mt_rate

# observed number of segregating sites
obs_central <- pan_data$SegregatingSites[1]
obs_western <- pan_data$SegregatingSites[2]

# number of simulations
nsim <- 1000

# tolerance (closest 10% of simulations)
tol <- 0.1 

# sample size (sample size for each subspecies is 5)
n <- ind_num # or is it 10?

# define the total mutation rate for the locus under study
mutrate <- mut_rate  # mutation rate per site
nsites <- num_sites # number of segregating sites
nloci <- 1 # number of loci



### 2. Define the prior distribution & parameter values

# Define the prior distribution for effective population size (Ne) and sample randomly from it
prior<- runif(nsim, 10, 100000) #param



### 3. Simulate data and compute summary statistics

# this is the same for both subspecies, because they have the same sample size, number of sites and mutation rate

# Simulate the segregating sites
simS <- numeric(nsim) # empty array to save the simulated segregating sites for each simulation 

# Call sim.tree.mut function in a loop to simulate data and compute summary statistics
for(i in 1:nsim) {
  # Simulate data 
  muttree <- sim.tree.mut(sample=n, current=prior[i], ancestral = prior[i], time = 0, nrep = nloci, mu=mutrate, L=nsites)
  
  # Get the number of segregating sites from simulated data
  simS[i] <- ncol(muttree$seg_sites[[1]])

}



### 4. Calculate the distance between the simulated data and the observed data.

dist_sumstat_central <- abs(simS-obs_S_central)
dist_sumstat_western <- abs(simS-obs_S_western)


### 5. Reject the parameter values that result in large distances than the tolerance distance

tol_dst_central <- quantile(dist_sumstat_central, tol) # get the tolerance distance for the central subspecies
tol_dst_western <- quantile(dist_sumstat_western, tol) # get the tolerance distance for the western subspecies



### 6. Retain the parameter values that result in small distances.

closest_sims_central <- which(dist_sumstat_central<tol_dst_central) # get the closest simulations for the central subspecies
closest_sims_western <- which(dist_sumstat_western<tol_dst_western) # get the closest simulations for the western subspecies



### 7.  Visualize the results

# Plot the joint distribution of the prior and summary statistics 
abc_plot_rejection(simS, prior, dist_sumstat_central, tol) 
abc_plot_rejection(simS, prior, dist_sumstat_western, tol) 

# plot the posterior distribution of Ne
hist(prior[closest_sims_central], breaks=20, main="Central", xlab="Ne")
hist(prior[closest_sims_western], breaks=20, main="Western", xlab="Ne")

# Get the summary of the posterior distribution
summary(prior[closest_sims_central])
summary(prior[closest_sims_western])

