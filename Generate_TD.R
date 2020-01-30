### Generating the transformation-based dendrogram
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

library(igraph)
library(picante)

# Set sample name
Sample_Name = "Dataset_Name"

# Load in data
setwd("/path/to/ICR_data") # Set directory containing bulk transformation files
peak.2.peak = read.csv(list.files(pattern = "*peak.2.peak.csv")) # Load in the peak.2.peak file for all/bulk transformations
num.trans = read.csv(list.files(pattern = "*num.peak.trans.csv")) # Load in the num.peak.trans file for all/bulk transformations

# Altering data to match igraph structure
# peak.2.peak file
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak')] = 'id'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.x')] = 'from'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.y')] = 'to'
colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'Trans.name')] = 'type'

# num.trans file
colnames(num.trans)[which(colnames(num.trans) == 'peak')] = 'id'


peak.2.peak = peak.2.peak[,-which(colnames(peak.2.peak) %in% c('Dist','Dist.plus','Dist.minus'))]

peak.2.peak$weight = 1

head(peak.2.peak)

print(length(which(!peak.2.peak$from %in% num.trans$id)))
print(length(which(!peak.2.peak$to %in% num.trans$id)))
print(date())

peak.2.peak = peak.2.peak[,c('from','to','type','weight','sample')]
num.trans = num.trans[,c('id','num.trans.involved.in','sample')]

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

# Generate the UPGMA tree - a neighbor joining one will not work with this large of an object
tree = as.phylo(hclust(as.dist(net.dist), method = "average"))

# Store tree
write.tree(tree, paste(Sample_Name, "_TD_UPGMA.tre", sep = ""))
