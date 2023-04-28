# AAB 2023
# Individual Assignment: FST and Permutations
# Author: Cátia Antunes
# Date: 25-04-2023

##########################################################################################################################
source("utilfunctions.r")
library(hierfstat)

##########################################################################################################################

# 1. Function that implements a permutation test to assess if the FST
# between two populations from which the samples were taken is different
# from zero. This function takes as arguments two genotype matrices, and
# the number of permutations. The output is the p value.

# function that performs a permutation test 
permutation_test <- function(matrix1, matrix2, n_permutations) {
  # create a matrix of genotypes by combining the two matrices of genotypes so that the columns represent the SNPs and the rows represent the individuals
  genotypes <- cbind(matrix1, matrix2)
  
  # calculate the observed FST value for the populations using the getFst() function from the utilfunctions.r script
  fst_observed <- getFst(matrix1, matrix2)
  
  # create an empty vector with the length of the number of permutations to store the FST values from the permutations
  fst_permutations <- numeric(n_permutations)
  
  # loop through the specified number of permutations
  for (i in 1:n_permutations) {
    
    # randomly permute the indices of the individuals using the sample() function so that the individuals are randomly assigned to populations
    permuted_indices <- sample(1:ncol(genotypes))
    
    # split the permuted indices into two groups of indices, one for each population 
    permuted_indices1 <- permuted_indices[1:(ncol(matrix1))] # the first half of the permuted indices
    permuted_indices2 <- permuted_indices[(ncol(matrix1)+1):ncol(genotypes)] # the second half of the permuted indices
    
    # create permuted matrices of genotypes by selecting the columns of the original matrices using the permuted indices
    permuted_matrix1 <- matrix1[, permuted_indices1] # ???
    permuted_matrix2 <- matrix2[, permuted_indices2] # ???
    
    # combine the permuted matrices of genotypes into a single matrix so that the columns represent the SNPs and the rows represent the individuals
    #permuted_genotypes <- cbind(permuted_matrix1, permuted_matrix2)
    
    # calculate the FST for the permuted genotypes using the fst() function from the hierfstat package 
    fst_permutations[i] <- getFst(permuted_matrix1, permuted_matrix2)
  }
  
  # calculate the p-value as the proportion of permuted FST values that are greater than the observed FST
  p_value <- sum(fst_permutations >= fst_observed) / n_permutations
  
  return(p_value)
}

##########################################################################################################################

# 2. Applying the function defined above to assess the significance of the pairwise FST values
# between all pairwise comparisons for the dataset of Henn et al. (2015), using 1000 permutations, 
# and a significance level of 0.05.

# import hierfstat package
library(hierfstat)

# read the genotype matrix
genotype <- as.matrix(read.table("./AAB_Project/Exercise_1/Hennetal_genotypeMatrix.geno", header=T, stringsAsFactors = F, na.strings = "NA"))

# Get the number of individuals
nr_ind <- ncol(genotype)
nr_ind #54

# Get the number of SNPs
nr_snps <- nrow(genotype)
nr_snps # 10000

# Get the label of all individuals
ind <- colnames(genotype)
ind

# Read the file that contains the information regarding individuals and corresponding population
ind_pop <- read.table("./AAB_Project/Exercise_1/IndPopID.txt", header=T, stringsAsFactors = F, na.strings = "NA")
ind_pop

# Get the population labels
pop_names <- unique(ind_pop$popID)
pop_names
length(pop_names) #7

# Get a vector with the index of each individual according to the population
pop_indx <- numeric(nrow(ind_pop)) # initialize variable
for(individual in 1:length(ind)) { # loop through each individual in geno matrix
  # get the index of the line in indpop matrix
  # that corresponds to the name of ind
  row_individual <- which(ind_pop$IndID==ind[individual])
  # get the population of individual ind
  pop_of_ind <- ind_pop$popID[row_individual]
  # get the index of the population 
  pop_indx[individual] <- which(pop_names==pop_of_ind)
}
pop_indx # validated pop_indx[20] gives the population of individual 20

# Get a list with the index of individuals grouped by population
# using for loops
sample_indx <- list() 
for(population in 1:length(pop_names)) { # loop through pop.names
  # apply the function to detect which individuals belong to each pop
  sample_indx[[population]] <- which(pop_indx==population)
}
sample_indx #validated sample_indx[7] gives all the individuals from population 7
sample_indx[1]


#######

# 3. Implementation

# create a list to store the genotype matrices for each population
populations <- list()

# Split the genotype matrix into population genotype matrices, 1 matrix per population based on sample index
for (i in 1:length(pop_names)) {
  populations[[i]] <- genotype[, sample_indx[[i]]]
}
length(populations) #7

populations[1]


# # Compute the pairwise FST
# pairfst <- matrix(NA, ncol=length(pop_names), nrow=length(pop_names))
# # go through each pair of populations
# for(i in 1:(length(pop_names)-1)) {
#   for(j in (i+1):length(pop_names)) {
#     # call the function to compute the pairwise FST between each pair of populations
#     pairfst[i,j] <-  getFst(genotype[,sample_indx[[i]]],genotype[,sample_indx[[j]]])       
#   }
# }

library(hierfstat)


# define all pairwise comparisons
comparisons <- combn(length(populations), 2)

# loop through pairwise comparisons and calculate p-value
for (i in seq_len(ncol(comparisons))) {
  matrix1 <- populations[[comparisons[1,i]]]
  matrix2 <- populations[[comparisons[2,i]]]
  p_value <- permutation_test(matrix1, matrix2, n_permutations=1000)
  if (p_value < 0.05) {
    print(paste0("Pairwise FST between populations ", comparisons[1,i], " and ", comparisons[2,i], " is significant (p-value = ", round(p_value, 4), ")"))
  } else {
    print(paste0("Pairwise FST between populations ", comparisons[1,i], " and ", comparisons[2,i], " is not significant (p-value = ", round(p_value, 4), ")"))
  }
}

# keep getting error Error in matrix1[, permuted_indices1] : subscript out of bounds