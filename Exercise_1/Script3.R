#
source("utilfunctions.r")
library(fsthet)

# read the genotype matrix
geno <- as.matrix(read.table("./Hennetal_genotypeMatrix.geno", header=T, stringsAsFactors = F, na.strings = "NA"))
# check if we read genotype matrix correctly
str(geno) # matrix with nsites rows and nind columns
# How many individuals? How many sites? Individuals =54 ; sites = 10000

# get the label of individuals for each column of geno
ind <- colnames(geno)

# read a file with information about individuals and corresponding population
indpop <- read.table("./IndPopID.txt", header=TRUE, stringsAsFactors = FALSE)

# Population names according to distance from San - vector
pop.names <- c("San","Mbuti","Mozabite","Pathan","Cambodian","Yakut","Maya")
pop.names

# Get a vector with the index of population for each individual
# using for loops
popinds <- numeric(nrow(indpop)) # initialize variable
for(i in 1:length(ind)) { # loop through each individual in geno matrix
  # get the index of the line in indpop matrix
  # that corresponds to the name of ind
  row_i <- which(indpop$IndID==ind[i])
  # get the population of individual ind
  pop_of_ind <- indpop$popID[row_i]
  # get the index of the population 
  popinds[i] <- which(pop.names==pop_of_ind)
}
popinds

  
# Get a list with the index of individuals that belong to each population
# using for loops
sample_inds <- list() 
for(i in 1:length(pop.names)) { # loop through pop.names
  # apply the function to detect which individuals belong to each pop
  sample_inds[[i]] <- which(popinds==i)
}



##########
permutation_test_fst <- function(population1, population2, n_permutations) {
  # Combine the two populations
  combined_data <- cbind(population1, population2)
  
  # Calculate the observed FST
  n_snps <- nrow(population1)
  n_ind <- ncol(combined_data)
  obs_fst <- calc_fst(population1, population2)
  
  # Initialize a vector to store the FST values from the permutations
  perm_fst <- numeric(n_permutations)
  
  # Perform the permutations
  for (i in 1:n_permutations) {
    # Shuffle the individuals between the two populations
    shuffle <- sample(n_ind)
    shuffled_data <- combined_data[, shuffle]
    
    # Split the shuffled data back into two populations
    shuffled_population1 <- shuffled_data[, 1:ncol(population1)]
    shuffled_population2 <- shuffled_data[, (ncol(population1)+1):ncol(shuffled_data)]
    
    # Calculate the FST for the shuffled populations
    perm_fst[i] <- calc_fst(shuffled_population1, shuffled_population2)
  }
  
  # Calculate the p-value as the proportion of permuted FST values that are
  # greater than or equal to the observed FST
  p_value <- sum(perm_fst >= obs_fst) / n_permutations
  
  return(p_value)
}

# Helper function to calculate FST
calc_fst <- function(population1, population2) {
  n_snps <- nrow(population1)
  n_ind1 <- ncol(population1)
  n_ind2 <- ncol(population2)
  
  p_bar <- (colSums(population1) + colSums(population2)) / (n_ind1 + n_ind2)
  
  num <- sum((colMeans(population1) - colMeans(population2))^2) - (1/n_snps)*sum(p_bar*(1-p_bar)*((n_ind1-n_ind2)^2)/(n_ind1+n_ind2))
  den <- sum(p_bar*(1-p_bar)) - (1/n_snps)*sum(p_bar*(1-p_bar)*((n_ind1-n_ind2)^2)/(n_ind1+n_ind2))
  
  fst <- num / den
  
  return(fst)
}
