# AAB 2023
# Individual Assignment: Introduction to Approximate Bayesian Computation (ABC)
# Author: Cátia Antunes
# Date: 25-04-2023

##########################################################################################################################

# Load the packages
library("locfit")
library("scrm")
library("ape")
library(scrm)
source("coalfunctions.r")

##########################################################################################################################

setwd("./Exercise_2") # set working directory

# Load data
pan_data <- read.table("SegregatingSitesChimp.txt", header = TRUE)

pop <- pan_data$Population # load population
pop
ind_num <- pan_data$NumIndividuals # load number of individuals
ind_num
seg_sites <- pan_data$SegregatingSites # load segregating sites
seg_sites
num_sites <- pan_data$NumSites # load numbber of sites
num_sites
mt_rate <- pan_data$MtRatePerSite # load mutation rate per site
mt_rate

# observed sumstat
obs_S <- c(126, 64) # number of segregating sites for the Central and Western subspecies

# number of simulations
nsim <- 1000

# Define sample size
n <- 5 # number of individuals sampled from each population
nsites <- 50000 # number of sites in the genome
nloci <- 1 # number of loci used to estimate the mutation rate

# Define the total mutation rate for the locus under study
mutrate <- 1.2e-8  

# Define the prior distribution for the parameter (effective population size Ne)
prior <- function(n) runif(n, 10, 100000) 

# Define the summary statistic function to compute the difference between observed and simulated number of segregating sites
sumstat <- function(data) abs(data - obs_S)

# Define the distance function as the Euclidean distance between the observed and simulated summary statistics
distance <- function(s, s.sim) sqrt(sum((s - s.sim)^2))

# Simulate data and perform ABC
abc.res <- abc(nsim = nsim, prior = prior, sumstat = sumstat, distance = distance,
               param.names = "Ne", n = n, nsites = nsites, nloci = nloci, mutrate = mutrate)

# get the posterior distribution of the effective population size for each subspecies
posterior_Central <- abc.res$posterior[,"Ne"][abc.res$posterior[,"param"]=="Central"]
posterior_Western <- abc.res$posterior[,"Ne"][abc.res$posterior[,"param"]=="Western"]
