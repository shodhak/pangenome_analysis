rm(list = ls())

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
#Assign a new variable "origin" that tells whether the sample is from Hungary, Israel, or CDC
sample_info$origin <- "CDC"
sample_info$origin[which(substr(sample_info$yale_code,1,1) == "H")] <- "Hungary"
sample_info$origin[which(substr(sample_info$yale_code,1,1) == "D")] <- "Israel"

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

#Remove 0 from serotype names in tmp_sample_info
tmp_sample_info$pneumocat_serotype <- as.character(tmp_sample_info$pneumocat_serotype)
#Remove the 0 when it is at the starting at the name of serotype
for(i in 1:length(tmp_sample_info$pneumocat_serotype)){
  if (substr(tmp_sample_info$pneumocat_serotype[i],1,1) == "0"){
    tmp_sample_info$pneumocat_serotype[i] = substring(tmp_sample_info$pneumocat_serotype[i],2)
  }
  else {
    tmp_sample_info$pneumocat_serotype[i] = tmp_sample_info$pneumocat_serotype[i]
  }
}
tmp_sample_info$pneumocat_serotype <- as.factor(tmp_sample_info$pneumocat_serotype)

#REGRESSION ANALYSIS
#The gene_presence_absence.csv file contains the gene names, their sizes, number fo isolates that contain them,
#and also their annotation. Extract that information into a data frame
gene_info <- read.csv("gene_presence_absence.csv")
gene_info <- gene_info[,1:14]

#Load pangenome_gene_presence_absence.csv 
#This dataframe contains sample wise gene presence absence which was obtained as output from Roary
pangene_pa <- read.csv("pangenome_gene_presence_absence.csv")
#Make gene names as row names
rownames(pangene_pa) <- pangene_pa$Gene
#Delete first column
pangene_pa <- pangene_pa[,-1]

#Transform the gene presence absence data into a new data frame
t_pangene_pa <- data.frame(t(pangene_pa))

#Extract the samples for which the serotypes are not "NT", and are verified with PneumoCAT
t_pangene_pa <- t_pangene_pa[which(rownames(t_pangene_pa) %in% tmp_sample_info$sample_name),]
#Load polysaccharide composition data
ps_composition <- read.csv("ps_composition.csv")
#Add serotype column to pangenome data and make it a character variable
t_pangene_pa$serotype <- as.character(tmp_sample_info$pneumocat_serotype[match(rownames(t_pangene_pa),tmp_sample_info$sample_name)])
#Move serotype to first column
t_pangene_pa <- t_pangene_pa[,c(ncol(t_pangene_pa), 1:(ncol(t_pangene_pa)-1))]
#Add GlcA info to pangenome data
t_pangene_pa$GlcA <- ps_composition$GlcA[match(t_pangene_pa$serotype, ps_composition$Serotype)]
#Move GlcA to second column
t_pangene_pa <- t_pangene_pa[,c(1, ncol(t_pangene_pa), 2:(ncol(t_pangene_pa)-1))]

#The information aboue presence of GlcA is not available for some of the serotypes
GlcA_NA <- t_pangene_pa$serotype[which(is.na(t_pangene_pa$GlcA) == TRUE)]
#Drop the serotypes for which GlcA is NA
t_pangene_pa <- t_pangene_pa[-which(t_pangene_pa$serotype %in% GlcA_NA),]

#Run regression analysis for GlcA with respect to the genes
#Following is the code that was previously used for regression of acetyl and alcohol data. Modify it for current analysis.

#The following if_factor variable counts the number of levels for each column in a dataframe
if_factor <- sapply(fdat2, function(x) nlevels(as.factor(x)))
#To test the number of columns with more than one levels
length(which(if_factor>1))
#To test the number of columns with only one level
length(which(if_factor==1))

#Logistic Regression for Acetyl
reg.ace.nb1 <- vector("list", ncol(fdat2)) 
aic.ace.nb1 <- 1: ncol(fdat2)
pred.ace.nb1<-matrix(nrow = nrow(fdat2), ncol = ncol(fdat2) )
for (i in 1:ncol(fdat2)){
  if(nlevels(as.factor(fdat2[,i])) > 1){
    reg.ace.nb1[[i]]<-glm(as.factor(st_access1$acetyl) ~ as.factor(fdat2[,i]), family = binomial)
    aic.ace.nb1[[i]]<-reg.ace.nb1[[i]]$aic
    pred.ace.nb1[,i]<-predict(reg.ace.nb1[[i]])
  }
  else{
    #AIC values cannot be calculated since there is only one level
    reg.ace.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    aic.ace.nb1[[i]]<-1000
    #Predicted value will be same as existing value
    pred.ace.nb1[,i]<-fdat2[1,i]
  }
}
#Compute weights for each position
wgt_acetyl <- exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1)))/sum(exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1))))

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
wgt1_acetyl <- wgt_acetyl
wgt1_acetyl[wgt1_acetyl == wgt1_acetyl[1]] <- 0
pos1 <- colnames(fdat2)
#Make a dataframe containing positions and weights
plot_data_ace <- data.frame(pos1,wgt1_acetyl) 
names(plot_data_ace) <- c("Position","Weight")
plot_data_ace$Product <- ps_gene$Product[ps_gene$Position %in% plot_data_ace$Position]

#Write csv for plot_data to make a shinyapp for the plot (Refer to app.R)
write.csv(plot_data_ace,"plot_data_acetyl.csv")

#Plot weights with respect to position
plot_data_ace$Position <- factor(plot_data_ace$Position, as.character(plot_data_ace$Position))
ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight), group=1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")

#Magnify the effects
#^0.005brings weight values closer to 1
#^10 magnifies the effect
#Plotly makes interactive plot
plot_acetyl <- ggplot(data=plot_data_ace, aes(x=plot_data_ace$Position, y=(plot_data_ace$Weight^0.005)^10, group=1, text=plot_data_ace$Product)) +
  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
  scale_x_discrete(breaks = plot_data_ace$Position[which(plot_data_ace$Weight!=0)]) +
  labs(x = "Position", y = "Weight (adjusted)")
ggplotly(plot_acetyl, tooltip = c("x","text"))

#One of the positions has disproportainately high value of weight
#plot_data$Position[which(plot_data$Weight==max(plot_data$Weight))] #1023758

#Remove the weight for that position
#plot_data1 <- plot_data[-which(plot_data$Position==1023758),]

#Replot
#plot_data1$Position <- factor(plot_data1$Position, as.character(plot_data1$Position))
#ggplot(data=plot_data1, aes(x=plot_data1$Position, y=(plot_data1$Weight^0.005)^10, group=1)) +
#  geom_line() +theme(axis.text.x = element_text(angle=90, hjust=1, size = 5),panel.background = element_blank()) +
#  scale_x_discrete(breaks = plot_data1$Position[which(plot_data1$Weight!=0)]) +
#  labs(x = "Position", y = "Weight (adjusted)")


