### A script to analyze null modeling results

# Load in required libraries
require(reshape2); require(ggplot2); require(ggthemes)
require(vegan)


#### Load in data ####

setwd("/path/to/null_model_results")

# Load in bNTI (either active or inactive bNTI go here)
char = read.csv("Shoreline_Active40_MCD_bNTI_999.csv", row.names = 1)
trans = read.csv("Shoreline_Active40_TD_bNTI_999.csv", row.names = 1)
weig = read.csv("Shoreline_Active40_TWCD_bNTI_999.csv", row.names = 1)


# ######################## #
#### Data Preprocessing ####
# ######################## #

# Load in poor calibrations
poor.cal = read.csv("Shoreline_Poorly_Calibrated_Samples.csv")

# Removing poorly calibrated samples
char = char[-which(row.names(char) %in% poor.cal$samples), -which(colnames(char) %in% poor.cal$samples)]
trans = trans[-which(row.names(trans) %in% poor.cal$samples), -which(colnames(trans) %in% poor.cal$samples)]
weig = weig[-which(row.names(weig) %in% poor.cal$samples), -which(colnames(weig) %in% poor.cal$samples)]

rm("poor.cal")

# Reflecting null matrices
char[upper.tri(char)] = t(char)[upper.tri(char)]
trans[upper.tri(trans)] = t(trans)[upper.tri(trans)]
weig[upper.tri(weig)] = t(weig)[upper.tri(weig)]


# ############################## #
#### ggplot objects for later ####
# ############################## #

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

color.list = c("#4575b4", "#e2e2ac", "#d73027")


# ####################### #
#### Testing out means ####
# ####################### #
# All samples
char.mean = data.frame(bNTI = apply(char, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Characteristics", stringsAsFactors = F)
tran.mean = data.frame(bNTI = apply(trans, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Transformations", stringsAsFactors = F)
weig.mean = data.frame(bNTI = apply(weig, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Weighted", stringsAsFactors = F)

merge.mean = rbind(char.mean, tran.mean, weig.mean)
merge.mean$Type[grep("River", row.names(merge.mean))] = "Surface Water"

all.mean.stats = data.frame(Dend = unique(merge.mean$Dend), W = NA, p.value = NA, stringsAsFactors = F)

for(i in 1:length(all.mean.stats$Dend)){
  temp.stats = merge.mean[which(merge.mean$Dend %in% all.mean.stats$Dend[i]),]
  temp.mwu = wilcox.test(bNTI~Type, data = temp.stats)
  
  all.mean.stats$W[i] = temp.mwu$statistic
  all.mean.stats$p.value[i] = temp.mwu$p.value
}

ggplot(data = merge.mean, aes(x = Type, y = bNTI))+
  geom_boxplot(aes(group = Type, color = Type))+
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  ggtitle("All Samples Means")+
  facet_grid(Dend~., scales = "free_y")+
  scale_color_stata()+
  hori_x_theme

# Within-group samples
char.temp = char
char.temp[grep("River_at", row.names(char.temp)), -grep("River_at", colnames(char.temp))] = NA
char.temp[-grep("River_at", row.names(char.temp)), grep("River_at", colnames(char.temp))] = NA

tran.temp = trans
tran.temp[grep("River_at", row.names(tran.temp)), -grep("River_at", colnames(tran.temp))] = NA
tran.temp[-grep("River_at", row.names(tran.temp)), grep("River_at", colnames(tran.temp))] = NA

weig.temp = weig
weig.temp[grep("River_at", row.names(weig.temp)), -grep("River_at", colnames(weig.temp))] = NA
weig.temp[-grep("River_at", row.names(weig.temp)), grep("River_at", colnames(weig.temp))] = NA

char.mean = data.frame(bNTI = apply(char.temp, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Characteristics", stringsAsFactors = F)
tran.mean = data.frame(bNTI = apply(tran.temp, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Transformations", stringsAsFactors = F)
weig.mean = data.frame(bNTI = apply(weig.temp, 2, mean, na.rm = T), Type = "Pore Water", Dend = "Weighted", stringsAsFactors = F)

merge.mean = rbind(char.mean, tran.mean, weig.mean)
merge.mean$Type[grep("River", row.names(merge.mean))] = "Surface Water"

within.mean.stats = data.frame(Dend = unique(merge.mean$Dend), W = NA, p.value = NA, stringsAsFactors = F)

for(i in 1:length(within.mean.stats$Dend)){
  temp.stats = merge.mean[which(merge.mean$Dend %in% within.mean.stats$Dend[i]),]
  temp.mwu = wilcox.test(bNTI~Type, data = temp.stats)
  
  within.mean.stats$W[i] = temp.mwu$statistic
  within.mean.stats$p.value[i] = temp.mwu$p.value
}

ggplot(data = merge.mean, aes(x = Type, y = bNTI))+
  geom_boxplot(aes(group = Type, color = Type))+
  geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
  ggtitle("Within Groups Means")+
  facet_grid(Dend~., scales = "free_y")+
  scale_color_stata()+
  hori_x_theme

rm("merge.mean", "temp.mwu", "temp.stats", "char.mean", 
   "tran.mean", "weig.mean", "char.temp", "tran.temp", "weig.temp")


# ##################################### #
#### Preprocessing Data for plotting ####
# ##################################### #

# Melting data
char = melt(as.matrix(char)); char$Tree = "Characteristics"
trans = melt(as.matrix(trans)); trans$Tree = "Transformations"
weig = melt(as.matrix(weig)); weig$Tree = "Weighted"

# Removing null values
char = char[!is.na(char$value),]
trans = trans[!is.na(trans$value),]
weig = weig[!is.na(weig$value),]

# Determining ecological processes fractionation
eco.proc = data.frame(hom.sel = rep(NA, 3), stoch = NA, var.sel = NA, 
                      row.names = c("Characteristics", "Transformation", "Weighted"))
eco.proc$hom.sel = c((length(which(char$value <= -2))/length(char$value)), 
                     (length(which(trans$value <= -2))/length(trans$value)), 
                     (length(which(weig$value <= -2))/length(weig$value)))
eco.proc$stoch = c((length(which(abs(char$value) < 2))/length(char$value)), 
                   (length(which(abs(trans$value) < 2))/length(trans$value)), 
                   (length(which(abs(weig$value) < 2))/length(weig$value)))
eco.proc$var.sel = c((length(which(char$value >= 2))/length(char$value)), 
                     (length(which(trans$value >= 2))/length(trans$value)), 
                     (length(which(weig$value >= 2))/length(weig$value)))

eco.proc = melt(as.matrix(eco.proc))
eco.proc$value = (eco.proc$value*100)

ggplot(data = eco.proc, aes(x = "", y = value, fill = Var2))+
  geom_bar(stat = "identity")+
  coord_polar("y", start=0)+
  ggtitle("Ecological Processes")+
  xlab(NULL) + ylab(NULL)+
  facet_grid(Var1~.)+
  scale_fill_manual(values = c("dodgerblue4", "goldenrod3", "firebrick4"))+
  theme_bw()+
  theme(axis.text.x=element_blank())
  

# # Combining data
# null.data = rbind(char, trans, weig)
# null.data = null.data[!is.na(null.data$value),]
# 
# # Renaming samples to simpler names
# sample_names = as.character(null.data$Var1)
# sample_names = gsub("Stegen_", "", sample_names); sample_names = gsub("_Nov.*$", "", sample_names)
# null.data$Var1 = sample_names
# 
# sample_names = as.character(null.data$Var2)
# sample_names = gsub("Stegen_", "", sample_names); sample_names = gsub("_Nov.*$", "", sample_names)
# null.data$Var2 = sample_names
# 
# rm("char", "trans", "weig", "sample_names")
# 
# # Adding Location Information
# null.data$Type = "Pushpoint"
# null.data$Type[grep("River_at", null.data$Var2)] = "River"
# 
# 
# # ###################### #
# #### Plotting results ####
# # ###################### #
# 
# ### All Samples
# # Density plots
# ggplot(data = null.data, aes(x = value, fill = Tree))+
#   geom_density(alpha = 0.6)+
#   geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("All Samples")+
#   scale_fill_manual(values = color.list)+
#   hori_x_theme
# 
# # Boxplots by sample
# ggplot(data = null.data, aes(x = Var2, y = value))+
#   geom_boxplot(aes(group = Var2, color = Type))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("All Samples")+
#   facet_grid(Tree~., scales = "free_y")+
#   scale_color_stata()+
#   vert_x_theme
# 
# # Boxplots by group
# ggplot(data = null.data, aes(x = Type, y = value))+
#   geom_boxplot(aes(group = Type, color = Type))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("All Samples")+
#   facet_grid(Tree~., scales = "free_y")+
#   scale_color_stata()+
#   hori_x_theme
# 
# # All sample stats
# all.stats = data.frame(Dendrogram = unique(null.data$Tree), F.stat = NA, aov.P = NA, MWU.W = NA, MWU.P = NA)
# 
# for(i in 1:length(all.stats$Dendrogram)){
#   temp.stats = null.data[which(null.data$Tree %in% all.stats$Dendrogram[i]),]
#   
#   temp.ano = summary(aov(value~Type, data = temp.stats))
#   temp.mwu = wilcox.test(value~Type, data = temp.stats)
#   
#   all.stats$F.stat[i] = temp.ano[[1]]$`F value`[1]
#   all.stats$aov.P[i] = temp.ano[[1]]$`Pr(>F)`[1]
#   all.stats$MWU.W[i] = temp.mwu$statistic
#   all.stats$MWU.P[i] = temp.mwu$p.value
#   
# }
# 
# 
# ### Pushpoint only
# temp = null.data[-c(grep("River_at", null.data$Var1), grep("River_at", null.data$Var2)),]
# 
# # Density plots
# ggplot(data = temp, aes(x = value, fill = Tree))+
#   geom_density(alpha = 0.3)+
#   geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("Pushpoint Only")+
#   scale_fill_manual(values = color.list)+
#   hori_x_theme
# 
# # Boxplots by sample
# ggplot(data = temp, aes(x = Var2, y = value))+
#   geom_boxplot(aes(group = Var2))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("Pushpoint Only")+
#   facet_grid(Tree~., scales = "free_y")+
#   vert_x_theme
# 
# 
# ### River only
# temp = null.data[intersect(grep("River_at", null.data$Var1), grep("River_at", null.data$Var2)),]
# 
# # Density plots
# ggplot(data = temp, aes(x = value, fill = Tree))+
#   geom_density(alpha = 0.3)+
#   geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("River Only")+
#   scale_fill_manual(values = color.list)+
#   hori_x_theme
# 
# # Boxplots by sample
# ggplot(data = temp, aes(x = Var2, y = value))+
#   geom_boxplot(aes(group = Var2))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("River Only")+
#   facet_grid(Tree~., scales = "free_y")+
#   vert_x_theme
# 
# 
# ### Without Cross-comparisons
# temp = null.data[-c(grep("River_at", null.data$Var1), grep("River_at", null.data$Var2)),]
# temp = rbind(temp, null.data[intersect(grep("River_at", null.data$Var1), grep("River_at", null.data$Var2)),])
# 
# # Density plots
# ggplot(data = temp, aes(x = value, fill = Tree))+
#   geom_density(alpha = 0.3)+
#   geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("Without Cross-comparisons")+
#   scale_fill_manual(values = color.list)+
#   hori_x_theme
# 
# # Boxplots by sample
# ggplot(data = temp, aes(x = Var2, y = value))+
#   geom_boxplot(aes(group = Var2, color = Type))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   xlab(label = NULL)+
#   ggtitle("Without Cross-comparisons")+
#   facet_grid(Tree~., scales = "free_y")+
#   scale_color_stata()+
#   vert_x_theme
# 
# # Boxplots by group
# ggplot(data = temp, aes(x = Type, y = value))+
#   geom_boxplot(aes(group = Type, color = Type))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("Without Cross-comparisons")+
#   facet_grid(Tree~., scales = "free_y")+
#   scale_color_stata()+
#   hori_x_theme
# 
# # Stats
# within.stats = data.frame(Dendrogram = unique(temp$Tree), F.stat = NA, aov.P = NA, MWU.W = NA, MWU.P = NA)
# 
# for(i in 1:length(within.stats$Dendrogram)){
#   temp.stats = temp[which(temp$Tree %in% within.stats$Dendrogram[i]),]
#   
#   temp.ano = summary(aov(value~Type, data = temp.stats))
#   temp.mwu = wilcox.test(value~Type, data = temp.stats)
#   
#   within.stats$F.stat[i] = temp.ano[[1]]$`F value`[1]
#   within.stats$aov.P[i] = temp.ano[[1]]$`Pr(>F)`[1]
#   within.stats$MWU.W[i] = temp.mwu$statistic
#   within.stats$MWU.P[i] = temp.mwu$p.value
#   
# }
# 
# 
# ### Cross-comparisons only
# temp = null.data[intersect(grep("^ICR", null.data$Var1), grep("River_at", null.data$Var2)),]
# temp = rbind(temp, null.data[intersect(grep("River_at", null.data$Var1), grep("^ICR", null.data$Var2)),])
# 
# # Density plots
# ggplot(data = temp, aes(x = value, fill = Tree))+
#   geom_density(alpha = 0.3)+
#   geom_vline(xintercept = c(-2,2), color = "red", lty = 2)+
#   ggtitle("Cross-Comparison Only")+
#   scale_fill_manual(values = color.list)+
#   hori_x_theme
# 
# # Boxplots
# ggplot(data = temp, aes(x = Var2, y = value))+
#   geom_boxplot(aes(group = Var2, color = Type))+
#   geom_hline(yintercept = c(-2,2), color = "red", lty = 2)+
#   xlab(label = NULL)+
#   ggtitle("Cross-Comparison Only")+
#   facet_grid(Tree~., scales = "free_y")+
#   scale_color_stata()+
#   vert_x_theme
