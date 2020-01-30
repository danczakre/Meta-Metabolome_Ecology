### Calculates RCBC for ICR data
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013

match.tree = T # This ensures that the data matches the tree used for bNTI
tree.type = "Tree_Type"

library(vegan); library(picante)

rc.reps = 9999 # Set the number of null iterations

# ############################### #
#### Data Loading and cleaning ####
# ############################### #

setwd("/path/to/ICR_data") # Working directory for FT-ICR data
data = read.csv("Processed_Data.csv", row.names = 1) # Importing the organismal data  
tree = read.tree("Metabolite_Dendrogram.tre") # Importing the tree

# Converting to presence/absence
if(max(data) > 1){
  data[data > 1] = 1
}

# Matching data to tree, if desired
if(match.tree){
  phylo = match.phylo.data(tree, data)
  
  data = phylo$data
  tree = phylo$phy
  
  rm("phylo")
}

# ######################## #
#### Running Raup-Crick ####
# ######################## #

# Setting up the data to mesh with the Stegen et al. code
spXsite = t(data)

# Count number of sites and total species richness across all plots (gamma)
n_sites = nrow(spXsite)
gamma = ncol(spXsite)

# Build a site by site matrix for the results, with the names of the sites in the row and col names:
results = matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))

# Make the spXsite matrix into a new, pres/abs. matrix:
spXsite.inc = ceiling(spXsite/max(spXsite))

# Create an occurrence vector- used to give more weight to widely distributed species in the null model
occur = apply(spXsite.inc, MARGIN=2, FUN=sum)

# Create an abundance vector- used to give more weight to abundant species in the second step of the null model
abundance = apply(spXsite, MARGIN=2, FUN=sum)

# Loops through every pairwise comparison, generating null results

for(null.one in 1:(nrow(spXsite)-1)){
  for(null.two in (null.one+1):nrow(spXsite)){
    
    null_bray_curtis<-NULL
    
    for(i in 1:rc.reps){
      
      # Generates two empty communities of size gamma
      com1<-rep(0,gamma)
      com2<-rep(0,gamma)
      
      # Add observed number of species to com1, weighting by species occurrence frequencies
      com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
      
      # Again for com2
      com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
      
      null.spXsite = rbind(com1,com2); # Null.spXsite
      
      # Calculates the null Bray-Curtis
      null_bray_curtis[i] = vegdist(null.spXsite, method='bray');
      
    }; # End of the null loop
    
    # Unlisting the parallel list
    null_bray_curtis = unlist(null_bray_curtis)
    
    # Calculates the observed Bray-Curtis
    obs.bray = vegdist(spXsite[c(null.one,null.two),], method='bray');
    
    # How many null observations is the observed value tied with?
    num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
    
    # How many null values are smaller than the observed *dissimilarity*?
    num_less_than_in_null = sum(null_bray_curtis<obs.bray);
    
    rc = ((num_less_than_in_null + (num_exact_matching_in_null)/2)/rc.reps) # This variation of rc splits ties
    
    rc = (rc-.5)*2 # Adjusts the range of the  Raup-Crick caclulation to -1 to 1
    
    results[null.two,null.one] = round(rc,digits=2); # Stores rc into the results matrix
    
    print(c(null.one,null.two,date())); # Keeps track of position
    
  }; # End of inner loop
  
}; # End of outer loop

rc.results = as.dist(results) # Converts results into a distance matrix

write.csv(as.matrix(rc.results), paste0("Shoreline_", tree.type, "_RCBC_", rc.reps, ".csv"), quote = F)

rm('spXsite.inc')