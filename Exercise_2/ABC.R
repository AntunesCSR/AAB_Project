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

getwd()
setwd("C:/Users/antun/OneDrive/Documents/Education/MRes/Year 1 - Semester 2/AAB/Final Project/AAB_Project/Exercise_2") # set working directory

##########################################################################################################################


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
mt_rate <- 1.2e-8 # had to manually change the value of mutation rate because it was not being read correctly
mt_rate

# observed number of segregating sites
obs_S_central <- pan_data$SegregatingSites[1]
obs_S_western <- pan_data$SegregatingSites[2]

# number of simulations
nsim <- 1000

# tolerance (closest 10% of simulations)
tol <- 0.1 
tol05 <- 0.05
tol40 <- 0.4

# sample size for each population
n <- 10

# define the total mutation rate for the locus under study
mutrate <- mt_rate  # mutation rate per site
nsites <- num_sites # number of segregating sites
nloci <- 1 # number of loci



### 2. Define the prior distribution & parameter values

# Define the prior distribution for effective population size (Ne) and sample randomly from it
param <- runif(nsim, 10, 100000) 



### 3. Simulate data and compute summary statistics

# this is the same for both subspecies, because they have the same sample size, number of sites and mutation rate

# Simulate the segregating sites
simS <- numeric(nsim) # empty array to save the simulated segregating sites for each simulation 

# Call sim.tree.mut function in a loop to simulate data and compute summary statistics
for(i in 1:nsim) {
  # Simulate data 
  muttree <- sim.tree.mut(sample=n, current=param[i], ancestral = param[i], time = 0, nrep = nloci, mu=mutrate, L=nsites)
  
  # Get the number of segregating sites from simulated data
  simS[i] <- ncol(muttree$seg_sites[[1]])

}



### 4. Calculate the distance between the simulated data and the observed data.

dist_sumstat_central <- abs(simS-obs_S_central)
dist_sumstat_western <- abs(simS-obs_S_western)



### 5. Reject the parameter values that result in large distances than the tolerance distance

tol_dst_central <- quantile(dist_sumstat_central, tol) # get the tolerance distance for the central subspecies (tolerance 10%)
tol_dst_western <- quantile(dist_sumstat_western, tol) # get the tolerance distance for the western subspecies  (tolerance 10%)

tol05_dst_central <- quantile(dist_sumstat_central, tol05) # get the tolerance distance for the central subspecies (tolerance 5%)
tol05_dst_western <- quantile(dist_sumstat_western, tol05) # get the tolerance distance for the western subspecies  (tolerance 5%)

tol40_dst_central <- quantile(dist_sumstat_central, tol40) # get the tolerance distance for the central subspecies (tolerance 40%)
tol40_dst_western <- quantile(dist_sumstat_western, tol40) # get the tolerance distance for the western subspecies  (tolerance 40%)



### 6. Retain the parameter values that result in small distances.

closest_sims_central <- which(dist_sumstat_central<tol_dst_central) # get the closest simulations for the central subspecies (tolerance 10%)
closest_sims_western <- which(dist_sumstat_western<tol_dst_western) # get the closest simulations for the western subspecies (tolerance 10%)

closest_sims_central_05 <- which(dist_sumstat_central<tol05_dst_central) # get the closest simulations for the central subspecies (tolerance 5%)
closest_sims_western_05 <- which(dist_sumstat_western<tol05_dst_western) # get the closest simulations for the western subspecies (tolerance 5%)

closest_sims_central_40 <- which(dist_sumstat_central<tol40_dst_central) # get the closest simulations for the central subspecies (tolerance 40%)
closest_sims_western_40 <- which(dist_sumstat_western<tol40_dst_western) # get the closest simulations for the western subspecies (tolerance 40%)



### 7.  Visualize the results

# Plot and save the joint distribution of the prior and summary statistics for both subspecies
png("central_joint_distribution_tol10.png", width=600, height=800)
obs_S <- obs_S_central 
abc_plot_rejection(simS, param, dist_sumstat_central, tol) # (tolerance 10%)
dev.off()

png("central_joint_distribution_tol05.png", width=600, height=800)
obs_S <- obs_S_central 
abc_plot_rejection(simS, param, dist_sumstat_central, tol05) # (tolerance 5%)
dev.off()

png("central_joint_distribution_tol40.png", width=600, height=800)
obs_S <- obs_S_central
abc_plot_rejection(simS, param, dist_sumstat_central, tol40) # (tolerance 40%)
dev.off()

png("western_joint_distribution_tol10.png", width=600, height=800)
obs_S <- obs_S_western
abc_plot_rejection(simS, param, dist_sumstat_western, tol) # (tolerance 10%)
dev.off()

png("western_joint_distribution_tol05.png", width=600, height=800)
obs_S <- obs_S_western
abc_plot_rejection(simS, param, dist_sumstat_western, tol05) # (tolerance 5%)
dev.off()

png("western_joint_distribution_tol40.png", width=600, height=800)
obs_S <- obs_S_western
abc_plot_rejection(simS, param, dist_sumstat_western, tol40) # (tolerance 40%)
dev.off()


# Plot and save the posterior distribution for both subspecies
png("central_posterior_histograms.png", width=600, height=800)
hist(param[closest_sims_central], breaks=20, main="Central", xlab="Ne")
dev.off()

png("western_posterior_histograms.png", width=600, height=800)
hist(param[closest_sims_western], breaks=20, main="Western", xlab="Ne")
dev.off()


# Get the summary of the posterior distribution for both subspecies and save to file

summary(param[closest_sims_central]) # (tolerance 10%)
param_central <- as.array(summary(param[closest_sims_central]))
write.csv(param_central, file = "central_posterior.csv")

summary(param[closest_sims_central_05]) # (tolerance 5%)
param_central <- as.array(summary(param[closest_sims_central]))
write.csv(param_central, file = "central_posterior.csv")

summary(param[closest_sims_central_40]) # (tolerance 40%)
param_central <- as.array(summary(param[closest_sims_central]))
write.csv(param_central, file = "central_posterior.csv")

summary(param[closest_sims_western]) # (tolerance 10%)
param_western <- as.array(summary(param[closest_sims_western]))
write.csv(param_western, file = "western_posterior.csv")

summary(param[closest_sims_western_05]) # (tolerance 5%)
param_western <- as.array(summary(param[closest_sims_western]))
write.csv(param_western, file = "western_posterior.csv")

summary(param[closest_sims_western_40]) # (tolerance 40%)
param_western <- as.array(summary(param[closest_sims_western]))
write.csv(param_western, file = "western_posterior.csv")
