### Calculating null bMNTD variants
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013


range = 1:999
Sample_Name = "Dataset_Name"
tree_type = "Tree_Type"

#-----------------#

library(vegan)
library(picante)

###################################
#### Data Loading and cleaning ####
###################################

setwd("/path/to/ICR_data") # Working directory for FT-ICR data
data = read.csv("Processed_Data.csv", row.names = 1) # Importing the organismal data  
tree = read.tree("Metabolite_Dendrogram.tre") # Importing the tree

# Creating necessary directories

if(!dir.exists(paste0(tree.type, "_Null_Results"))){
  dir.create(paste0(tree.type, "_Null_Results"))
}


####################################
#### Beginning the bNTI Process ####
####################################

# Converting to presence/absence
data[data>1] = 1

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree, data)

# Running cophenetic outside of the for-loop
coph = cophenetic(phylo$phy)

# Calculating the bMNTD for 999 random distributions
print(paste(date(), " - Start for loop"))

for(i in range){
  bMNTD.rand = as.matrix(comdistnt(t(phylo$data), taxaShuffle(coph), abundance.weighted = F, exclude.conspecifics = F))
  write.csv(bMNTD.rand, paste(tree_type, "_Null_Results/FTICR_", Sample_Name, "_", tree_type, "_bMNTD_rep", i, ".csv", sep = ""), quote = F)
  rm("bMNTD.rand")
  
  print(c(date(),i))
} # Performing the calculations on using the OTU table but with randomized taxonomic affiliations

print(paste(date(), " - End for loop"))
