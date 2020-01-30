### Generating the molecular characteristics dendrogram
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

options(digits = 10)

require(phangorn) # For tree based functions
require(ggtree) # For tree visualization
require(vegan) # For vegdist


# ################## #
#### Load in data ####
# ################## #

# Set sample name
Sample_Name = "Dataset_Name"

# Load in data
setwd("/path/to/ICR_data")
mol = read.csv(list.files(pattern = "*Mol.csv"), row.names = 1) # Load in molecular data

### Ensuring that isotopic peaks are removed
if(length(which(colnames(data) %in% "C13")) > 0){
  w = row.names(data)[which(data$C13 > 0)]
  
  if(length(w) > 0){
    data = data[-which(row.names(data) %in% w),]
  }
  
  rm("w")
}

# Removing peaks that have no formula assignments
mol = mol[-which(mol$MolForm %in% NA),]

# Setting objects for useful parameters
Mol.Info = mol[,c("C", "H", "O", "N", "S", "P", "DBE", "AI_Mod", "kdefect"), drop = F]
Mol.Ratio = mol[,c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")]


# ##################### #
#### Generating tree ####
# ##################### #

# Pairwise distance between peaks
Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info)) # Generating a distance matrix based upon the provided parameters

# Create tree
tree = as.phylo(hclust(vegdist(Mol.Info, "euclidean"), "average")) # Converting the distance matrix into a tree

# Quick visualization (with Van Krevelen compound classes)
mol = cbind(row.names(mol), mol) # Peak names need to be the first column for ggtree to work

col = colorRampPalette(c("dodgerblue4", "goldenrod3", "firebrick3"))(length(unique(mol$bs1_class)))
ggtree(tree, layout = "circular") %<+% mol + geom_tippoint(aes(color = bs1_class), na.rm = T) +
  scale_color_manual(values = col)+
  theme(legend.position = "right")

# Writing tre
write.tree(tree, paste(Sample_Name, "_MCD_UPGMA.tre", sep = ""))
