### Generating the transformation-weighted characteristics dendrogram
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

options(digits = 10)

require(phangorn) # For tree based functions
require(ggtree) # For tree visualization
require(vegan) # For vegdist
library(igraph) # For network functionality

# ################## #
#### Load in data ####
# ################## #

# Set sample name
Sample_Name = "Dataset_Name"

# Load in data
setwd("/path/to/ICR_data")
mol = read.csv(list.files(pattern = "*Mol.csv"), row.names = 1) # Load in the molecular data
peak.2.peak = read.csv(list.files(pattern = "*peak.2.peak.csv")) # Load in the peak.2.peak file for all/bulk transformations
num.trans = read.csv(list.files(pattern = "*num.peak.trans.csv")) # Load in the num.peak.trans file for all/bulk transformations

# Before doing anything, I'm removing peaks that have an isotope signature
if(length(which(colnames(mol) %in% "C13")) > 0){
  w = row.names(mol)[which(mol$C13 > 0)]
  
  if(length(w) > 0){
    mol = mol[-which(row.names(mol) %in% w),]
  }
  
  rm("w")
}

# Removing peaks that have no formula assignments
mol = mol[-which(mol$MolForm %in% NA),]

# Setting objects for useful parameters
Mol.Info = mol[,c("C", "H", "N", "P", "O", "S", "DBE", "AI_Mod", "kdefect"), drop = F]
Mol.Rat = mol[,c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")]

# ####################### #
#### Cleaning the data ####
# ####################### #

## Altering transformation to match igraph structure
# peak.2.peak file
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak')] = 'id'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.x')] = 'from'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.y')] = 'to'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'Trans.name')] = 'type'
peak.2.peak$weight = 1

# num.trans file
colnames(num.trans)[which(colnames(num.trans) == 'peak')] = 'id'

# Ensuring the two data sets match
print(length(which(!peak.2.peak$from %in% num.trans$id)))
print(length(which(!peak.2.peak$to %in% num.trans$id)))

# Reordering the data
peak.2.peak = peak.2.peak[,c('from','to','type','weight','sample')]
num.trans = num.trans[,c('id','num.trans.involved.in','sample')]


# ########################### #
#### Determining distances ####
# ########################### #

### Pairwise molecular distance between peaks
Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info)) # Generating a distance matrix based upon the provided parameters
mol.dist = as.matrix(vegdist(Mol.Info, "euclidean"))

### Determining transformation distance
# Creating the network
net = graph_from_data_frame(d=peak.2.peak, vertices=num.trans, directed=F)
rm("peak.2.peak", "num.trans")

# The distances command is much better than the similarity measurement
net.dist = distances(net)

# Finding clusters and determining the distance in the largest
clus = clusters(net)
max.clus = which(clus$csize %in% max(clus$csize)) # Finding the largest cluster
max.clus = names(clus$membership)[which(clus$membership %in% max.clus)] # Finding the members of the largest cluster

net.dist = net.dist[max.clus, max.clus] # Setting the net dist to that size only; only lost ~4000 peaks by doing this with the dereplicated HJ-Andrews set

# Need to normalize the dissimiarlity to 0-1
net.dist = (net.dist-min(net.dist))/(max(net.dist)-min(net.dist))

### Parsing down the data
q = which(row.names(net.dist) %in% row.names(mol.dist)) # Net dist is generate in the Merged_EdgeNode_Files script - I probably will incorporate it here (or something)
net.dist.data = net.dist[q, q] # Matching the network distance to the molecular information

q = which(row.names(mol.dist) %in% row.names(net.dist.data))
mol.dist.data = mol.dist[q, q] # Matching the molecular information to the network distance

# Weighting and tree generation
weighted.dist = mol.dist.data*net.dist.data # Weighting the tree


# ######################### #
#### Generating the tree ####
# ######################### #

# Creating tree
tree = as.phylo(hclust(as.dist(weighted.dist), method = "average"))

# Quick visualization (with Van Krevelen compound classes)
mol = cbind(row.names(mol), mol) # Peak names need to be the first column for ggtree to work
ggtree(tree, layout = "circular") %<+% mol + geom_tippoint(aes(color = bs1_class), na.rm = T) +
  theme(legend.position = "right")

# Writing the tree
write.tree(tree, paste(Sample_Name, "_Weighted_All-Trans_UPGMA.tre", sep = ""))
