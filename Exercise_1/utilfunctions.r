# Definition of functions .....................................................................

# IMAGE.SCALE
# This function creates a color scale for use with the image() function. 
# Input parameters should be consistent with those used in the corresponding image plot. 
# The "horiz" argument defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

#' EXPHET_SITE
#' computes the expected heterozygosity for a given site
#' INPUT:
#'  @param geno_site vector of size nind with the genotypes coded as 0,1,2 and NA (missing data)
#' OUTPUT:
#'   expected heterozygosity for a given site
ExpHet_site <- function(geno_site) {
  
  # get the number of individuals with data
  ngenecopies <- 2*sum(!is.na(geno_site))
  
  # get the frequency of the alternative allele
  freq <- sum(geno_site, na.rm=T)/ngenecopies
  
  # output the expected heterozygosity
  he <- (ngenecopies/(ngenecopies-1))*2*freq*(1-freq)
  he
}

#' EXPHET
#' computes the average expected heterozygosity across all sites
#' INPUT:
#'  @param geno_matrix  matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#' OUTPUT:
#'  @return average expected heterozygosity
ExpHet <- function(geno_matrix) {
  # compute the heterozygosity for each site
  exp_het <- apply(geno_matrix, 1, function(row) {ExpHet_site(row)}) 
  # compute the mean across sites
  mean(exp_het)
}

#' EXPHET_MEAN
#' 'computes the average expected heterozygosity across all sites
#' INPUT:
#'  @param geno_site  matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#'  @param sample_inds list where each entry is a vector of size nind_pop with the index of individuals belonging to a given pop
#'  @param pop_names vector of size npop with the names of each population
#' OUTPUT:
#'  @retunr average expected heterozygosity
ExpHet_mean <- function(geno_matrix, sample_inds, pop.names) {
  
  # initialize matrix to save heterozygosity for each SNP and for each pop
  # het is a matrix with nsnps rows and npop columns
  exp_het <- matrix(NA,nrow=nrow(geno_matrix), ncol=length(pop.names))
  # go through each population
  for(i in 1:length(pop.names)) {
    # compute the heterozygosity by looking at a subset of individuals from geno matrix
    exp_het[,i] <- apply(geno_matrix[,sample_inds[[i]]], 1, function(row) {ExpHet_site(row)}) 
  }
  # check the output
  #str(exp_het)
  
  # get the mean het for each pop
  mean_exp_het <- colMeans(exp_het, na.rm=TRUE)
  names(mean_exp_het) <- pop.names
  mean_exp_het
}

#' FST
#' compute FST according to Hudson's estimator following Bathia
#' INPUT:
#'  @param geno1  matrix with nsites x nind1 with the genotype data for population 1
#'  @param geno2  matrix with nsites x nind2 with the genotype data for population 2
#' OUTPUT
#'  @return FST between the two populations
getFst <- function(geno1, geno2) {
  
  # compute the sample size for each site
  ss1 <- 2*rowSums(!is.na(geno1))
  ss2 <- 2*rowSums(!is.na(geno2))
  
  # compute allele frequency for each site
  freq1 <- rowSums(geno1, na.rm=T)/ss1
  freq2 <- rowSums(geno2, na.rm=T)/ss2
  
  # compute the terms p1(1-p1) and p2(1-p2)
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  
  
  # compute the square of the difference among the allele frequenies
  pdiffsquare <- (freq1-freq2)^2
  
  # compute the numerator
  numerator <- pdiffsquare-(p1/(ss1-1))-(p2/(ss2-1))
  
  # compute the denominator
  denominator <- (freq1*(1-freq2))+(freq2*(1-freq1))
  
  
  # output FST estimators
  sum(numerator, na.rm=T)/sum(denominator, na.rm=T)
}





#' PLOTFST
#' function to plot the pairwise FST values
#' INPUT:
#'  @param fst matrix with npop rows and npop columns with the pairwise FST
#'  @param pnames vector of size npop with the population labels
#' OUTPUT
#'  @return plot with the pairwise FST
plotFst <- function(fst, pnames) {
  # layout defines that the plot area has two regions with different sizes
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(0.8,0.2))
  # use Rcolor brewer to define a color scale
  nclasses <- colorRampPalette(brewer.pal(8,"Oranges"))(15)
  # define the breakpoints
  breaksplot <- seq(0, max(fst, na.rm=T)*1.02,length.out = length(nclasses)+1)
  # plot the image
  image(1:nrow(fst), 1:nrow(fst), t(fst), breaks=breaksplot, xlab="", ylab="", col=nclasses, main="Pairwise FST", xaxt="n", yaxt="n")
  axis(1, at=c(1:nrow(fst)), labels=pnames)
  axis(2, at=c(1:nrow(fst)), labels=pnames)
  image.scale(t(fst), zlim=range(breaksplot), col = nclasses, breaks=breaksplot, horiz=FALSE, ylab="Pairwise FST", cex.lab=1.2, cex.axis=1.2, cex.main=1.2)
}
