rm(list = ls())

#The gene_presence_absence.csv file contains the gene names, their sizes, number fo isolates that contain them,
#and also their annotation. Extract that information into a data frame
gene_info <- read.csv("gene_presence_absence.csv")
gene_info <- gene_info[,1:14]

#Load pangenome_gene_presence_absence.csv 
#This dataframe contains sample wise gene presence absence
pangene_pa <- read.csv("pangenome_gene_presence_absence.csv")

#Load information about individual samples so that associations between the above datasets and the serotypes
#and MLST can be made. That will further help in performing logistic regression with respect to the 
#presence/absence of GlcA

#Load file resequencing_check.csv that contains sample names and numbers assigned in Adrienn's experiments
#These also contains output from PneumoCAT, and the column miscall represents whether the serotypes matched or not.
#Miscall value of 1 means a mismatch and value of 0 means matching.
#The mismatch were sent for resequencing so till that data come, these samples may be dropped from the analysis
sample_info <- read.csv("resequencing_check.csv")

#The above dataset DOES NOT contain information about the sample names Pool_DMW and Pool_hung that will be necessary to
#tie the results with results from sequencing data
#Load data for DMW and hun samples
dmw_samples <- read.csv("WGS_Isr_corrected.csv")
#Select first five columns, and the Name column that contains the Yale code for the samples
dmw_samples <- dmw_samples[,c(1:4,10)]
#Assign names to column
colnames(dmw_samples) <- c("miscall","corrected_well","assigned_serotype","pneumocat_serotype","yale_code")
#Add string "Pool_DMW_" to wells for dmw_samples
dmw_samples$sample_name <- paste("Pool_DMW_",dmw_samples$corrected_well,sep = "")
#In miscall column, change NAs to 0
dmw_samples$miscall[which(is.na(dmw_samples$miscall)==TRUE)] <- 0

#Prepare similar data frame for Pool_hun samples
hun_samples <- read.csv("hun_sequencing.csv")
#Assign same column names as for the Hungary dataset
colnames(hun_samples) <- c("miscall","corrected_well","assigned_serotype","pneumocat_serotype","yale_code")
#Add string "Pool_hung_" to wells for hun_samples
hun_samples$sample_name <- paste("Pool_hung_",hun_samples$corrected_well,sep = "")

#Combine back all data frames
sample_info <- rbind(dmw_samples, hun_samples)

#Drop the samples for which the assigned serotype does not match the one PneumoCAT (misscall value 1)
sample_info <- sample_info[-which(sample_info$miscall == 1),]
#Remove non-typable serotypes
sample_info <- sample_info[-which(sample_info$assigned_serotype == "NT"),]
sample_info <- sample_info[-which(sample_info$pneumocat_serotype == "Failed"),]

#Load panini gene presence absence data
panini_pa <- read.csv("panini-gene_presence_absence.Rtab.csv")
#Add serotype information to panini data
#Find the samples that are common between the two data frames
samples_common <- intersect(sample_info$sample_name, panini_pa$id)

#Extract data for common samples from above from both datasets
tmp_panini_pa <- panini_pa[panini_pa$id %in% samples_common,]
tmp_sample_info <- sample_info[sample_info$sample_name %in% samples_common,]
tmp_panini_pa$serotype <- tmp_sample_info$pneumocat_serotype[match(tmp_sample_info$sample_name, tmp_panini_pa$id)]

#Transfer data from tmp_panini_pa to a new data frame
panini_sero <- tmp_panini_pa

panini_sero$serotype <- as.character(panini_sero$serotype)
#Remove the 0 when it is at the starting at the name of serotype
for(i in 1:length(panini_sero$serotype)){
  if (substr(panini_sero$serotype[i],1,1) == "0"){
    panini_sero$serotype[i] = substring(panini_sero$serotype[i],2)
  }
  else {
    panini_sero$serotype[i] = panini_sero$serotype[i]
  }
}
panini_sero$serotype <- as.factor(panini_sero$serotype)
#Extract serogroup from serotyep and make a new variable
panini_sero$serogroup <- gsub("[A-Z]","",panini_sero$serotype)

#Save the panini_sero as csv for visualization in microreact
write.csv(panini_sero, "panini_sero.csv")

#load library ggplot2 to plot the panini output
library(ggplot2)
p <- ggplot(data = panini_sero) + geom_point(aes(x = x, y = y, color = serotype), size = 3, alpha = 0.9)
p

p <- ggplot(data = panini_sero) + geom_point(aes(x = x, y = y, color = serogroup), size = 3, alpha = 0.9)
p

