# AAB 2023
# Individual Assignment: FST and Permutations
# Author: Cátia Antunes
# Date: 25-04-2023

##########################################################################################################################

# Load packages

library(RColorBrewer)
source("utilfunctions.r")

getwd()
setwd("C:/Users/antun/OneDrive/Documents/Education/MRes/Year 1 - Semester 2/AAB/Final Project/AAB_Project/Exercise_1") # set working directory

##########################################################################################################################

# 1. Function that implements a permutation test to assess if the FST
# between two populations from which the samples were taken is different
# from zero. This function takes as arguments two genotype matrices, and
# the number of permutations. The output is the p value.

permutation_test <- function(geno1, geno2, nperm) {
  
  # combine the genotype matrices
  geno12 <- cbind(geno1, geno2)

  # get the number of individuals in each population
  n1 <- ncol(geno1)
  n2 <- ncol(geno2)

  # compute the observed FST
  obs_fst <- getFst(geno1, geno2)

  # create a vector to store the permuted FST values
  perm_fst <- numeric(nperm)
  
  # go through each permutation
  for (i in 1:nperm) {
    # permute the index of individuals between populations
    perm_inds <- sample.int(n1 + n2, n1 + n2, replace = FALSE)
    # get the permuted genotype matrices
    perm_geno1 <- geno12[, perm_inds[1:n1]] 
    perm_geno2 <- geno12[, perm_inds[(n1+1):(n1+n2)]]
    # compute the permuted FST
    perm_fst[i] <- getFst(perm_geno1, perm_geno2)
  }

  # compute the p-value
  pval <- sum(perm_fst >= obs_fst) / nperm
  return(pval)
}

##########################################################################################################################

# 2. Applying the function defined above to assess the significance of the pairwise FST values
# between all pairwise comparisons for the dataset of Henn et al. (2015), using 1000 permutations, 
# and a significance level of 0.05.

# read the genotype matrix
geno <- as.matrix(read.table("./Hennetal_genotypeMatrix.geno", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA"))

# read a file with information about individuals and corresponding population
indpop <- read.table("./IndPopID.txt", header = TRUE, stringsAsFactors = FALSE)

# Get a vector with the index of population for each individual
popnames <- unique(indpop$popID)
popnames
popinds <- match(indpop$popID, popnames)
popinds

# Get a list with the index of individuals that belong to each population
sample_inds <- split(1:ncol(geno), popinds)
sample_inds

# Initialize matrix to save pairwise FST
pairfst <- matrix(NA, nrow = length(popnames), ncol = length(popnames))

# go through each pair of populations and compute the pairwise FST
for (i in 1:(length(popnames) - 1)) {
  for (j in (i + 1):length(popnames)) {
    # call the function to compute the pairwise FST between each pair of populations
    pairfst[i, j] <- permutation_test(geno[, sample_inds[[i]]], geno[, sample_inds[[j]]], nperm = 1000)
    pairfst[j, i] <- pairfst[i, j]
  }
}

# Print and export the matrix with the pairwise FST values
colnames(pairfst) <- popnames
rownames(pairfst) <- popnames
pairfst
write.csv(pairfst, file = "./pairwiseFST.csv", quote = FALSE)

# Print the matrix with the pairwise p-values
colnames(pairpvalues) <- popnames
rownames(pairpvalues) <- popnames
pairpvalues

# Plot and save the pairwise FST values
png("Pairwise_FST.png", width = 800, height = 800)
plotFst(pairfst, popnames)
dev.off()
