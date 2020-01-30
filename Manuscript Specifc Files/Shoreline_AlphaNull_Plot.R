### Processing alpha diversity results

library(reshape2)
library(ggplot2); library(ggthemes)
library(ggpubr)

setwd("/path/to/AlphaDiv")


# ################################## #
#### Loading and rearranging data ####
# ################################## #

### Load in normal alpha diversity
trans = read.csv("Shoreline_TD_AlphaDiv.csv")
char = read.csv("Shoreline_MCD_AlphaDiv.csv")
weig = read.csv("Shoreline_TWCD_AlphaDiv.csv")

### Load in null model alpha diversity
trans.null = read.csv("Shoreline_TD_AlphaNull.csv")
char.null = read.csv("Shoreline_MCD_AlphaNull.csv")
weig.null = read.csv("Shoreline_TWCD_AlphaNull.csv")

# Adding in River-Pushpoint information
trans.null$Type = trans$Type; char.null$Type = char$Type; weig.null$Type = weig$Type


### Melting and then combining data data
# Melting
trans = melt(trans, id.vars = c("Sample", "Type")); trans$Tree = "Transformations"
char = melt(char, id.vars = c("Sample", "Type")); char$Tree = "Characteristics"
weig = melt(weig, id.vars = c("Sample", "Type")); weig$Tree = "Weighted"

trans.null = melt(trans.null, id.vars = c("Sample", "Type")); trans.null$Tree = "Transformations"
char.null = melt(char.null, id.vars = c("Sample", "Type")); char.null$Tree = "Characteristics"
weig.null = melt(weig.null, id.vars = c("Sample", "Type")); weig.null$Tree = "Weighted"

# Combining tree-based data
tree.data = rbind(trans, char, weig)
null.data = rbind(trans.null, char.null, weig.null)

rm("trans", "char", "weig", "trans.null", "char.null", "weig.null") # Cleaning up

# Converting ses.mpd/mntd to NRI/NTI respectively
null.data$value = null.data$value*-1


# ####################################### #
#### Plotting comparisons across trees ####
# ####################################### #

# Plotting "phylogenetically" informed metrics
for(var in unique(tree.data$variable)){
  w = which(tree.data$variable %in% var)
  
  print(
    ggplot(data = tree.data[w,], aes(x = Type, y = value, group = Tree))+
      geom_boxplot(aes(group = Type, color = Type))+
      facet_grid(Tree~., scales = "free_y")+
      scale_color_stata()+
      xlab(label = NULL)+
      ylab(label = var)+
      theme_bw()+
      theme(text = element_text(size = 14),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.border = element_rect(size = 1, colour = "black"),
            panel.grid = element_blank(),
            legend.position = "none")
  )
}

# Alpha stats
names = NULL

for(i in 1:length(unique(tree.data$Tree))){
  for(j in 1:length(unique(tree.data$variable))){
    names = c(names, paste(unique(tree.data$Tree)[i], unique(tree.data$variable)[j]))
  }
}

alpha.stats = data.frame(Group = names, W = NA, p.value = NA, stringsAsFactors = F)

for(i in 1:length(alpha.stats$Group)){
  temp = unlist(strsplit(alpha.stats$Group[i], split = " "))
  temp.stats = tree.data[which(tree.data$Tree %in% temp[1] & tree.data$variable %in% temp[2]),]
  temp.mwu = wilcox.test(value~Type, data = temp.stats)
  
  alpha.stats$W[i] = temp.mwu$statistic
  alpha.stats$p.value[i] = temp.mwu$p.value
}


rm("temp.stats", "temp.mwu", "temp")

# Plotting the null modeling results
for(var in unique(null.data$variable)){
  w = which(null.data$variable %in% var)
  
  print(
    ggplot(data = null.data[w,], aes(x = Type, y = value, group = Tree))+
      geom_boxplot(aes(group = Type, color = Type))+
      facet_grid(Tree~., scales = "free_y")+
      scale_color_stata()+
      xlab(label = NULL)+
      ylab(label = var)+
      theme_bw()+
      theme(text = element_text(size = 14),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            panel.border = element_rect(size = 1, colour = "black"),
            panel.grid = element_blank(),
            legend.position = "none")
  )
}

# Null stats
null.stats = data.frame(Group = c(paste(unique(null.data$Tree), unique(null.data$variable)[1]),
                                       paste(unique(null.data$Tree), unique(null.data$variable)[2])), 
                                       W = NA, p.value = NA, stringsAsFactors = F)

for(i in 1:length(null.stats$Group)){
  temp = unlist(strsplit(null.stats$Group[i], split = " "))
  temp.stats = null.data[which(null.data$Tree %in% temp[1] & null.data$variable %in% temp[2]),]
  temp.mwu = wilcox.test(value~Type, data = temp.stats)

  null.stats$W[i] = temp.mwu$statistic
  null.stats$p.value[i] = temp.mwu$p.value
}


rm("temp.stats", "temp.mwu", "temp")

### Plotting together
# Richness
w = which(tree.data$variable %in% "SR")
p1 =  ggplot(data = tree.data[w,], aes(x = Type, y = value, group = Tree))+
  geom_boxplot(aes(group = Type, color = Type))+
  facet_grid(Tree~., scales = "free_y")+
  scale_color_stata()+
  xlab(label = NULL)+
  ylab(label = "Richness")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

# Phylogenetic diversity
w = which(tree.data$variable %in% "PD")
p2 =  ggplot(data = tree.data[w,], aes(x = Type, y = value, group = Tree))+
  geom_boxplot(aes(group = Type, color = Type))+
  facet_grid(Tree~., scales = "free_y")+
  scale_color_stata()+
  xlab(label = NULL)+
  ylab(label = "DD")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

# MPD
w = which(tree.data$variable %in% "MPD")
p3 =  ggplot(data = tree.data[w,], aes(x = Type, y = value, group = Tree))+
  geom_boxplot(aes(group = Type, color = Type))+
  facet_grid(Tree~., scales = "free_y")+
  scale_color_stata()+
  xlab(label = NULL)+
  ylab(label = "MPD")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

# MNTD
w = which(tree.data$variable %in% "MNTD")
p4 =  ggplot(data = tree.data[w,], aes(x = Type, y = value, group = Tree))+
  geom_boxplot(aes(group = Type, color = Type))+
  facet_grid(Tree~., scales = "free_y")+
  scale_color_stata()+
  xlab(label = NULL)+
  ylab(label = "MNTD")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

# Combining the individual plots
p = ggarrange(p1, p2, p3, p4)
p
