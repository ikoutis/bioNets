# Here we do the TOM similarity experiment
#
# First choose the soft thresholding based on the scale free topology
# Calculate the adjacency matrix based 
#       (i.e. calculate the correlation matrix and then power it to the soft thresholding)
# Calculate the Topological Overlap Matrix (TOM)
# 
#
# Niloofar Aghaieabiane
# July 2019
#


#Loading required packages and libraries and setting the directory
setwd('Z:/Study/Research/Prof Zhi and Koutis/Similarity/Codes/Real Data - Similarity')
library(dynamicTreeCut)
library(moduleColor)
library(WGCNA)
options(stringsAsFactors = FALSE)



#########################################Loading the Expression and EntrezID Files##############################
#Print the GSE number
GSE <- readline(prompt = "Enter the GSE ")
print(paste("GSE",GSE, sep = ""))
#
#
#loading the expression without outlier file
#A sample by genes matrix
caption_title <- paste("Choose ", "GSE", GSE, " expression file", sep = "")
exp_NoOutlier_temp <- read.csv(choose.files(caption = caption_title), header = FALSE)
exp_NoOutlier_temp <- exp_NoOutlier_temp[-1,]
#Convert to numeric
exp_NoOutlier <- matrix(NaN, dim(exp_NoOutlier_temp)[1], dim(exp_NoOutlier_temp)[2])
for(i in  1:dim(exp_NoOutlier_temp)[1]){
  for(j in 1:dim(exp_NoOutlier_temp)[2]){
    exp_NoOutlier[i,j] <- as.numeric(exp_NoOutlier_temp[i,j])
  }
}


#Loading the gene EntrezID corresponding with expression file 
caption_title <- paste("Choose ", "GSE", GSE, " EntrezID file", sep = "")
EntrezID_NoOutlier_temp <- read.csv(choose.files(caption = caption_title), header = FALSE)
EntrezID_NoOutlier_temp <- EntrezID_NoOutlier_temp[-1,]
#Convert to numeric
EntrezID_NoOutlier <- -1
for(i in 1:length(EntrezID_NoOutlier_temp)){
  EntrezID_NoOutlier <- rbind(EntrezID_NoOutlier, as.numeric(EntrezID_NoOutlier_temp[i]))
}
EntrezID_NoOutlier <- EntrezID_NoOutlier[-1,]

#
#attach the gene EntrezID to the expression filecolumns
colnames(exp_NoOutlier) <- as.vector(EntrezID_NoOutlier)
#
#Data Set: a matrix of #samples by #genes

#########################################################################################################



####################################Constructing Adjacency Matrix########################################
#Here we calculate the adjacency matrix of the network based on TOM similarity
#
#########Prepration for soft thresholding so that the network become scale-free topology
#Choose a set of soft-thresholding powers
powers <- c(c(1:20))
seq(from = 12, to = 20, by = 2)

#Call the network topology analysis function
sft = pickSoftThreshold(exp_NoOutlier, networkType = "signed" ,powerVector = powers, verbose = 5)

#Plot the result
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

#Scale-free topolog fit index as a function of the soft-thresholding power
dev.new()
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft thresholding(power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = 'red')

#Finding the lowest soft threshodl such that the value of R^2 appraoaches to 0.9
abline(h = 0.90, col = "red")

#Enter the lowest threshold based on the plot
softpower <- readline(prompt = "Enter the lowest threshold for which the R^2 approaches to 0.9 ")
softpower <- as.numeric(softpower)

dev.new()
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = 'red')

#Calculating adjacency matrix
adj <- adjacency(exp_NoOutlier, type = "signed", power = softpower)

#Calculating tipological overlap matrix
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM
#returned dissTOM  a #genes*#genes matrix of dissimilarity based on TOM

####################################################################################################################




###################################################Clustering#######################################################
#Here we find the clusters in the network constructed based on TOM similarity

geneTree_TOM <- hclust(as.dist(dissTOM), method = "average")

#Plot the resulting clustering tree
dev.new()
sizeGrWindow(12,9)
plot(geneTree_TOM, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

#Finding the cluters
#we are oooking for large module
minModuleSize <- 30

#Module identifiacion using dynamic tree cut
dynamicMods_TOM <- cutreeDynamic(dendro = geneTree_TOM, method = "tree",
                                 minClusterSize = minModuleSize, deepSplit = FALSE)

#Returns the label of the cluster per gene
table(dynamicMods_TOM)

#Converting the numeric labels to colors
dynamicColors_TOM <- labels2colors(dynamicMods_TOM)
table(dynamicColors_TOM)

#Plot the dendogram and the color underneath
dev.new()
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_TOM, dynamicColors_TOM, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Dendro and Cluster Color based on TOM similarity")

####################################################################################################################



#############################################Gene Ontology Enrichment Analysis######################################
print("Gene Ontology Enrichment Starts...")

#Enter the annotation db for the data set
##annotation_db <- readline(prompt = "Enter the annotation library ")
annotation_db <- "mouse4302.db"

#Enter the number of best p-values by user
##num <- readline(prompt = "Enter the best number of p-values for ontology ")
##numPvalue <- as.numeric(num)
numPvalue <- 10

moduleLabels_TOM <- dynamicMods_TOM
#Finding the frequency of the genes per cluster
#It has two clumns:
#    #moduleLabels_STS: Numeric label of clusters
#    #Freq: Number of genes per cluster
clusters_label_freq_TOM <- as.data.frame(table(moduleLabels_TOM))

#Creating the sumary data frame for best numPvalue ontology outcomes of the clusters
clusters_GO_TOM <- data.frame(matrix(nrow = 0, ncol = 9))
summary_column_names <- c("clusterNum", "GOtype", "GOID" ,"Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term")
colnames(clusters_GO_TOM) <- summary_column_names

##########Calling gene ontology function per cluster
#Grouping the genes belong to the same cluster and then calling the gene ontology function
for(i in 1:dim(clusters_label_freq_TOM)[1]){
  #Extracting the i-th label
  label <- clusters_label_freq_TOM$moduleLabels_TOM[i]
  label <- droplevels(label)
  label <- as.numeric(levels(label))
  #
  ##label <- 11
  cluster <- -1
  for(j in 1:length(moduleLabels_TOM)){
    if(moduleLabels_TOM[j] == label){
      cluster = rbind(cluster, EntrezID_NoOutlier[j])
    }
  }
  cluster <- cluster[-1]
  if(length(cluster) != clusters_label_freq_TOM$Freq[i])
    stop("the number of the elements in cluster is not equal to the total number of the genes in the module")
  #cluster contain the gene EntrezIDs of the genes belonging to i-th cluster
  
  #Calling gene ontology function for summarizing the ontology of the cluster
  source('GOstats.R')
  
  alarm <- paste("GO for cluster", i, "starts...", sep = " " )
  print(alarm)
  res <- GO_analysis(EntrezID_NoOutlier, cluster, annotation_db, numPvalue, i, label )
  alarm <- paste(" GO for cluster", i, "is done...", sep = " " )
  print(alarm)
  
  rm(list = "cluster")
  
  
  
  #attaching the ontology of cluster i to the previous clusters ontologies
  #if("plyr" %in% rownames(installed.packages()) == FALSE) {install.package("plyr")}
  library(plyr)
  clusters_GO_TOM = rbind.fill(clusters_GO_TOM, res)
  
}

print("Gene Ontology Enrichment is done.")

######################################################################################################################




#############################################Saving The Results#######################################################
#Save the module colors and labels and the summary of annotation
caption_title <- paste("GSE", GSE, "_TOM_Experiment.RData" ,sep = "")
save(moduleLabels_TOM, geneTree_TOM, clusters_GO_TOM, clusters_label_freq_TOM, file = caption_title )
#Followiing items will be saved:
#     moduleLabels_TOM           =>  numeric label of the clusters
#     geneTree_TOM               =>  denrogram of the clustering
#     clusters_GO_TOM            =>  Summary of the gene ontology annotation
#     clusters_label_freq_TOM    => table of the clusters label and their corresponding frequncies
######################################################################################################################





