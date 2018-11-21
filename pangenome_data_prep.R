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

