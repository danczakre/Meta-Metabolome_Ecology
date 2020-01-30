### Detecting putative biochemical transformations across the entire dataset
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

library(dplyr)
library(tidyr)

options(digits=10) # Sig figs in mass resolution data

Sample_Name = "Dataset_Name"
resume = F # This switch controls whether or not you want to resume a previously started analysis


####################### 
### Loading in data ###
#######################

# Loading in ICR data
setwd("/path/to/ICR_data")
data = read.csv(list.files(path = ".", pattern = "*_Data.csv", full.names = T), row.names = 1) # Keeping data and mol-data seperate to ensure they are unaltered
mol = read.csv(list.files(path = ".", pattern = "*_Mol.csv", full.names = T), row.names = 1)

colnames(data) = paste("Sample_", colnames(data), sep="")

# Checking row names consistency
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect in your row names")
}

# Loading in transformations
trans.full =  read.csv("Transformation_Database.csv")
trans.full$Name = as.character(trans.full$Name)


###########################################
### Running through the transformations ###
###########################################

# pull out just the sample names
samples.to.process = colnames(data)

# error term
error.term = 0.000010

# Creating a presence matrix for all peaks observed within the dataset
bulk.peaks = as.data.frame(cbind(row.names(mol), 1))
colnames(bulk.peaks) = c("peak", "Sample_All")
bulk.peaks[,2] = as.numeric(bulk.peaks[,2])

# Grouping peaks by sample - in this case, this is just a formality to put data into the correct format for downstream analysis
Sample_Peak_Mat <- bulk.peaks %>% gather("sample", "value", -1) %>% filter(value > 0) %>% select(sample, peak)
Sample_Peak_Mat[,2] = as.numeric(as.character(Sample_Peak_Mat[,2]))
colnames(Sample_Peak_Mat) = c("sample", "peak.x")

# Loading in previously generated data (if it exists)
if(resume == T){
  
  # Listing exisiting files
  files = list.files(path = paste("/Users/danc783/Bulk Peak Transformations (do not backup)/", Sample_Name, sep = ""),
                     full.names = T)
  
  # Merging files
  peak.2.peak = NULL
  
  for(f in files){
    temp = read.csv(f)
    peak.2.peak = rbind(peak.2.peak, temp)
    rm("temp")
  }
  
  # Parsing down file name to peak count
  start_peak = gsub("/Users.*Distance_Results_","",files); start_peak = as.numeric(gsub("_for.*.csv","", start_peak)) 
  
  # Assigning start peak
  start_peak = max(start_peak)
  
  # Ensuring start peak and ID match
  id = files[grep(paste("_", start_peak, "_", sep = ""), files)]
  id = gsub("/Users.*peak_","",id); id = as.numeric(gsub(".csv","", id)) 
  
  if(as.numeric(id) == as.numeric(row.names(data)[start_peak])){
    print("Everything matches - moving on with the transformation analysis")
  } else {
    stop("Something went wrong with loading previous data - I recommend starting from the beginning.")
  }
  
  # Clean up
  rm("id", "files")
  
} else {
  peak.2.peak = NULL
  start_peak = 1
}

# Running a loop to compare each peak to each other peak
for(i in start_peak:(nrow(data)-1)){ # I cannot stress the importance of the "-1" here...
  
  # Creating a data matrix to ensur no repeat or negative differences
  Distance_Results = Sample_Peak_Mat[-1:-i,] # Removing all peaks up to, and including the current peak
  Distance_Results$peak.y = Sample_Peak_Mat$peak.x[i] # Setting the peak of interest
  Distance_Results$Dist = Distance_Results$peak.x - Distance_Results$peak.y # Finding the difference between all peaks and the peak of interest
  
  # Adding in error terms to the matrix
  Distance_Results$Dist.plus = Distance_Results$Dist + error.term
  Distance_Results$Dist.minus = Distance_Results$Dist - error.term
  Distance_Results$Trans.name = -999
  
  # Reorganizing the data to make it applicable with other scripts
  Distance_Results = Distance_Results[,c("sample", "Dist", "peak.x", "peak.y", "Dist.plus", "Dist.minus", "Trans.name")]
  
  for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
    
    mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
    if (length(mass.diff) > 1) { break() }
    Distance_Results$Trans.name[ which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)  ] = current.trans
  }
  
  # Removing differences that didn't match any transformation
  Distance_Results = Distance_Results[-which(Distance_Results$Trans.name == -999),]
  
  # Storing individual peak information if the script needs to stop
  if(!dir.exists(paste("/Users/danc783/Bulk Peak Transformations (do not backup)/", Sample_Name, sep = ""))){
    dir.create(paste("/Users/danc783/Bulk Peak Transformations (do not backup)/", Sample_Name, sep = ""))
  }
  
  write.csv(Distance_Results, paste0("/Users/danc783/Bulk Peak Transformations (do not backup)/", Sample_Name, "/Distance_Results_", i, "_for_peak_", row.names(data)[i], ".csv"),
            quote = F, row.names = F)
  
  # Building a larger peak.2.peak file
  peak.2.peak = rbind(peak.2.peak, Distance_Results)
  
  print(paste("Finished running through peak #", i, " on ", date(), sep = ""))
  
}

# Creating a num.trans file for network generation
peak.stack = as.data.frame(c(peak.2.peak$peak.x, peak.2.peak$peak.y)); head(peak.stack)
peak.profile = as.data.frame(tapply(X = peak.stack[,1], INDEX = peak.stack[,1], FUN = 'length' )); dim(peak.profile)
colnames(peak.profile) = 'num.trans.involved.in'
peak.profile$sample = Sample_Name
peak.profile$peak = row.names(peak.profile)

write.csv(peak.2.peak, paste(Sample_Name, "_All-Trans_peak.2.peak.csv", sep = ""), quote = F, row.names = F)
write.csv(peak.profile, paste(Sample_Name, "_All-Trans_num.peak.trans.csv", sep = ""), quote = F, row.names = F)
