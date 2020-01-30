### Script to analyze Shoreline data

tree.type = "Tree_Type" # Presets include "MCD", "TD", and "TWCD"
pull.outlier = T # This will remove poorly calibrated samples
alpha.div = F # This is a switch to calculate alpha diversity or not - can take a while
match.phy = F # Matches ICR data to those metabolites found on a given dendrogram

# ############################# #
#### Load required libraries ####
# ############################# #

require(vegan) # For broad ecology functions
require(reshape2); require(ggplot2); require(ggthemes) # For prettier graphs
require(picante); require(phytools); require(GUniFrac); require(pez) # Enables tree-based analyses


# ################## #
#### Load in data ####
# ################## #

setwd("/path/to/ICR_data")

# FREDA-processed ICR Data
data = read.csv("Processed_Shoreline_FTICR_Data_clean.csv", row.names = 1)
mol = read.csv("Processed_Shoreline_FTICR_Mol_clean.csv", row.names = 1)

# Metabolite tree (trees are found in the supplemental information of the manuscript)
if(tree.type == "MCD"){
  tree = read.tree("Molecular Characteristics Dendrogram - UPGMA.tre")
} 

if(tree.type == "TD"){
  tree = read.tree("Transformation-based Dendrogram - UPGMA.tre")
} 

if (tree.type == "TWCD"){
  tree = read.tree("Transformation-Weighted Characteristics Dendrogram - UPGMA.tre")
}


# Loading in poorly calibrated samples and removing them if necessary
if(pull.outlier == T){
  poor.cal = read.csv("Shoreline_Poorly_Calibrated_Samples.csv", stringsAsFactors = F)
  data = data[,-which(colnames(data) %in% poor.cal$samples)]
  
  if(length(which(rowSums(data) == 0)) > 0){
    mol = mol[-which(rowSums(data) == 0),]
    data = data[-which(rowSums(data) == 0),]
  }
}


# ############# #
#### Errors ####
# ############ #

# Checking row names consistency between molecular info and data
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect: the mol. info and peak counts don't match")
}

# Checking to ensure "FREDA_Processing.R" was run
if(length(which(mol$C13 == 1)) > 0){
  stop("Isotopic signatures weren't removed")
}

if(length(grep("QC_SRFAII", colnames(data))) > 0){
  stop("Suwannee River standards are still in the data")
}

if(max(data) > 1){
  print("Data is not presence/absence")
  data[data > 1] = 1
}


# ######################################### #
#### Data clean-up and factor processing #### 
# ######################################### #

# Shortening names
sample_names = colnames(data)
sample_names = gsub("Stegen_", "", sample_names); sample_names = gsub("_Nov.*$", "", sample_names)
colnames(data) = sample_names

rm("sample_names")

# Creating factors
meta = data.frame(Samples = colnames(data), Type = "Pushpoint", stringsAsFactors = F)
meta$Type[grep("River", meta$Sample)] = "River"

# Creating a sample order object
sample_order = c("ICR_1_pp1", "ICR_2_pp1", "ICR_3_pp1", "ICR_4_pp1", "ICR_5_pp1", "ICR_6_pp1", 
                 "ICR_6_pp2", "ICR_6_pp3", "ICR_7_pp1", "ICR_7_pp2", "ICR_7_pp3", "ICR_8_pp1", "ICR_8_pp2", "ICR_8_pp3",
                 "ICR_9_pp1", "ICR_9_pp2", "ICR_9_pp3", "ICR_10_pp1", "ICR_10_pp2", "ICR_10_pp3", "River_at_ICR2p1",
                 "River_at_ICR2p2", "River_at_ICR3", "River_at_ICR6", "River_at_ICR7", "River_at_ICR8", "River_at_ICR9",
                 "River_at_ICR10")

# Creating ggplot themes
hori_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

vert_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

# Rooting tree at midpoint (needed for UniFrac)
tree = midpoint.root(tree) # Tested other root locations, output seems similar


# ######################################### #
#### Analyzing peak properties by sample ####
# ######################################### #

# Creating boxplots for traits
prop = data.frame(Sample = NA, Type = NA, Property = NA, value = NA)

for(i in 1:ncol(data)){
  # Selecting on data for the current sample
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = mol[row.names(temp),]
  
  # Keeping interesting variables
  temp = data.frame(GFE = temp$GFE, AI_Mod = temp$AI_Mod)
  temp = temp[!is.na(temp$GFE),] # Removing NAs
  
  # Melting to get data into long format
  temp = melt(temp)
  temp = data.frame(Sample = colnames(data)[i], Type = meta$Type[i], Property = temp$variable, value = temp$value)
  
  # Adding it to the larger dataset
  prop = rbind(prop, temp)
  
  # Clean up
  rm("temp")
}

prop = prop[!is.na(prop$Sample),] # Removing linger NA row

# Forcing a specific factor order for samples
prop$Sample = factor(prop$Sample, levels = sample_order)

print(
  ggplot(data = prop, aes(x = Sample, y = value, group = Property))+
    geom_boxplot(aes(group = Sample, color = Type))+
    facet_grid(Property~., scales = "free_y")+
    scale_color_stata()+
    ggtitle("Molecular Properties")+
    vert_x_theme
)


# ######################### #
#### Ecoloigcal Analyses ####
# ######################### #

### Alpha diversity
if(alpha.div == T){
  
  # Matching data to selected tree (only necessary for alpha diversity as this is done for beta automatically)
  if(match.phy == T){
    phylo = match.phylo.data(tree, data)
    tree = phylo$phy
    data = phylo$data
    mol = mol[which(row.names(mol) %in% row.names(data)),]
    
    rm("phylo")
  }
  
  # Determinging diversity
  div = data.frame(Sample = colnames(data), Type = meta$Type, PD = NA, 
                   SR = NA, MPD = NA, MNTD = NA, VPD = NA, VNTD = NA) # Creaating empty data frame to store data
  
  div[,c("PD", "SR")] = pd(t(data), tree) # Faith's PD and species richness
  
  comp.comm = comparative.comm(tree, t(data)) # Creating a "pez" compatible object
  div[,c("MPD", "MNTD", "VPD", "VNTD")] = generic.metrics(comp.comm, metrics = c(.mpd, .mntd, .vpd, .vntd)) # Processing other phylogenetic a-diversity metrics 
  rm(comp.comm) # Removing this object as it is no longer necessary
  
  # Diversity significance
  div.sig = data.frame(MWU.stat = rep(NA, 6), p.value = NA, row.names = c("PD", "SR", "MPD",
                                                                          "MNTD", "VPD", "VNTD")) # Creating object to store signfiicance
  
  div.sig[1,] = c(wilcox.test(div$PD~div$Type)$statistic, wilcox.test(div$PD~div$Type)$p.value)
  div.sig[2,] = c(wilcox.test(div$SR~div$Type)$statistic, wilcox.test(div$SR~div$Type)$p.value)
  div.sig[3,] = c(wilcox.test(div$MPD~div$Type)$statistic, wilcox.test(div$MPD~div$Type)$p.value)
  div.sig[4,] = c(wilcox.test(div$MNTD~div$Type)$statistic, wilcox.test(div$MNTD~div$Type)$p.value)
  div.sig[5,] = c(wilcox.test(div$VPD~div$Type)$statistic, wilcox.test(div$VPD~div$Type)$p.value)
  div.sig[6,] = c(wilcox.test(div$VNTD~div$Type)$statistic, wilcox.test(div$VNTD~div$Type)$p.value)
  
  # Plotting diversity
  div = melt(div, id.vars = c("Sample", "Type")) # Converting to long format
  
  div$Sample = factor(div$Sample, levels = sample_order) # Reorganizing factors
  
  print(
    ggplot(data = div, aes(x = Sample, y = value, group = variable))+
      geom_bar(stat = "identity", aes(fill = Type))+
      facet_grid(variable~., scales = "free_y")+
      scale_fill_stata()+
      ggtitle("Alpha Diversity Metrics")+
      vert_x_theme
  ) # Plotting the data
} # Added this switch to save time while working with beta diversity metrics


### Beta diversity
# Creating stats object to store results
beta.stats = data.frame(Pseudo.F = rep(NA, 4), p.value = NA, row.names = c("Jaccard", "JacSub", "UniFrac", "bMNTD"))

# Creating distance matrix
dist = vegdist(x = t(data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)
beta.stats[1,] = c(adonis2(dist~meta$Type, permutations = 999)$`F`[1], adonis2(dist~meta$Type, permutations = 999)$`Pr(>F)`[1])

# Plotting Jaccard NMDS
nms = metaMDS(dist, try = 40) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores
nms$Type = meta$Type # Adding meta-data

print(
  ggplot(data = nms, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    ggtitle("Overall NMDS")+
    scale_color_stata()+
    hori_x_theme
) # Plotting NMS graph

### Phylogenetic beta-diversity
## Matching data to the provided tree
phylo = match.phylo.data(tree, data) # Matching ICR dataset to the tree

## NMDS on tree-subsetted data
# Creating distance matrix
dist = vegdist(x = t(phylo$data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)
beta.stats[2,] = c(adonis2(dist~meta$Type, permutations = 999)$`F`[1], adonis2(dist~meta$Type, permutations = 999)$`Pr(>F)`[1])

# Plotting Jaccard NMDS
nms = metaMDS(dist, try = 40) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores
nms$Type = meta$Type # Adding meta-data

print(
  ggplot(data = nms, aes(x = NMDS1, y = NMDS2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    scale_color_stata()+
    ggtitle(paste0("Subset NMDS - ", tree.type))+
    hori_x_theme
) # Plotting NMS graph

## Performing UniFrac analyses
uni = GUniFrac(otu.tab = t(phylo$data), tree = phylo$phy) # Calculating the UniFrac distance
uni = uni$unifracs[,,"d_UW"] # Only interested in the unweighted dataset
beta.stats[3,] = c(adonis2(as.dist(uni)~meta$Type, permutations = 999)$`F`[1], adonis2(as.dist(uni)~meta$Type, permutations = 999)$`Pr(>F)`[1])

uni.pcoa = ape::pcoa(uni) # Generating a principal coordinate analysis for the Uni distances
uni.scores = as.data.frame(uni.pcoa$vectors) # Getting scores
uni.scores$Type = meta$Type

print(
  ggplot(data = uni.scores, aes(x = Axis.1, y = Axis.2))+
    geom_point(aes(color = Type, shape = Type), size = 2.5)+
    xlab(label = paste0("PCoA1 (", round((uni.pcoa$values$Relative_eig[1]*100), 2), "%)"))+
    ylab(label = paste0("PCoA2 (", round((uni.pcoa$values$Relative_eig[2]*100), 2), "%)"))+
    scale_color_stata()+
    ggtitle(paste0("UniFrac - ", tree.type))+
    hori_x_theme
) # Plotting the UniFrac PCoA

# Writing out stats
write.csv(beta.stats, paste0("Shoreline_", tree.type, "_BetaDiv_Sig.csv"), quote = F)
