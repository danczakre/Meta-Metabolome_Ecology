### Relating bNTI metrics

library(vegan) # Needed for the Mantel test

# Load in data
icr = read.csv("/path/to/ICR_null_results/FTICR_bNTI_999.csv", row.names = 1)
bact = read.csv("/path/to/microbial_null_results/OTU_bNTI_999.csv", row.names = 1)


## Preprocessing data

# Cleaning ICR
poor.cal = read.csv("Shoreline_Poorly_Calibrated_Samples.csv")
icr = icr[-which(row.names(icr) %in% poor.cal$samples), -which(colnames(icr) %in% poor.cal$samples)]

# Cleaning microbial samples which can low sequence counts
remove =  c("ICR_6_pp2", "ICR_8_pp1", "ICR_8_pp2", "ICR_8_pp3", "ICR_7_pp1", "ICR_7_pp2", "ICR_9_pp2",
            "ICR_10_pp1", "ICR_10_pp2", "ICR_10_pp3", "River_at_ICR10", "River_at_ICR8")

w = which(row.names(bact) %in% remove)
bact = bact[-w, -w]

rm("remove", "w", "poor.cal")

# Reflecting null matrices
icr[upper.tri(icr)] = t(icr)[upper.tri(icr)]
bact[upper.tri(bact)] = t(bact)[upper.tri(bact)]

# Renaming samples in ICR data
colnames(icr) = gsub("Stegen_", "", colnames(icr)); colnames(icr) = gsub("_Nov.*$", "", colnames(icr))
row.names(icr) = gsub("Stegen_", "", row.names(icr)); row.names(icr) = gsub("_Nov.*$", "", row.names(icr))

row.names(icr) = gsub("ICR2p", "ICR2_", row.names(icr)); colnames(icr) = gsub("ICR2p", "ICR2_", colnames(icr))

# Matching datasets
w = which(row.names(icr) %in% row.names(bact))
icr = icr[w, w]

w = which(row.names(bact) %in% row.names(icr))
bact = bact[w, w]

icr = icr[order(row.names(icr)), order(colnames(icr))]
bact = bact[order(row.names(bact)), order(colnames(bact))]


## Correlating the null models
mant = mantel(bact, icr, method = "spearman", permutations = 9999)
mant

avg.icr = apply(icr, 2, function(x) mean(x, na.rm = T))
avg.bact = apply(bact, 2, function(x) mean(x, na.rm = T))

## Melting data for further testing
icr = melt(as.matrix(icr))
bact = melt(as.matrix(bact))

plot(bact$value, icr$value)
