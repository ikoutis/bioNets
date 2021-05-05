# Here we do the fifth experiment for the inner product similarity
#Experiment5: exper5
# Take sigmoid function on observation
# Take inner product on the observation
# Normalize the resulting inner product
# Apply scale free cretria
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
caption_title <- paste("Choose ", "GSE", GSE, "_expression file", sep = "")
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
caption_title <- paste("Choose ", "GSE", GSE, "_EntrezID file", sep = "")
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

rm(exp_NoOutlier_temp,EntrezID_NoOutlier_temp)
#######################################################################################################################


#############################################Constructing the adjacency matrix#########################################
#Here we construct the adjacency matrix for the experiment5

#Loading required functions
source('NeededFunctions_InnerProduct.R')

#Firstly, we take the sigmoid function on the observation
sig_exp <- sigOnObservation(exp_NoOutlier)
rm(exp_NoOutlier)

#Calculate the inner product of the t(sig_exp)%*%sig_exp
# The result is a matrix of #gene by #gene
transpose_sig_exp <- t(sig_exp)
innerProduct <- transpose_sig_exp %*% sig_exp
rm(transpose_sig_exp, sig_exp)
diag(innerProduct) <- 1
if(!isSymmetric(innerProduct))
  stop("The inner product is not symmetric")

#Normalized the inner product
norm_exp <- normalization(innerProduct)
rm(innerProduct)

####################################Applying Scale-Free Creteria######################################################
data <- norm_exp
powerVector <- c(seq(1, 10, by = 1), seq(12, 100, by = 2))
source('ScaleFree.R')
sft

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
     labels = powerVector, cex = cex1, col = 'red')

#Finding the lowest soft threshodl such that the value of R^2 appraoaches to 0.9
abline(h = 0.90, col = "red")

#Enter the lowest threshold based on the plot
softpower <- readline(prompt = "Enter the lowest threshold for which the R^2 approaches to 0.9 ")
softpower <- as.numeric(softpower)

dev.new()
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powerVector, cex = cex1, col = 'red')


#The adjacency matrix is created
adj_exper5 <- norm_exp^softpower
print("The adjacency matrix for experiment5 is created.")

#The dissimilarity of the adjacency
dissAdj_exper5 <- 1 - adj_exper5
######################################################################################################################



###########################################Clustering#################################################################
#Here we find the clusters based on the dissimilarity adjacency matrix
print("Start Clustering...")

geneTree_exper5 <- hclust(as.dist(dissAdj_exper5), method = "average")

#Plot the resulting clustering tree
dev.new()
sizeGrWindow(12,9)
plot(geneTree_exper5, xlab = "", sub = "", main = "Gene clustering for experiment5",
     labels = FALSE, hang = 0.04)

#Finding the cluters using dynamicTreeCut R-package

#we are oooking for large module
minModuleSize <- 30

#Use Dynamic Tree Cut algorithm
dynamicMods_exper5 <- cutreeDynamic(dendro = geneTree_exper5, method = "tree",
                                    minClusterSize = minModuleSize, deepSplit = FALSE)



#Returns the label of the cluster per gene
table(dynamicMods_exper5)

#Converting the numeric labels to colors
dynamicColors_exper5 <- labels2colors(dynamicMods_exper5)
table(dynamicColors_exper5)


#Plot the dendogram and the color underneath
dev.new()
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_exper5, dynamicColors_exper5, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Dendro and Cluster Color for experiment1")

#Clustering is done
print("Clustering for experiment5 is done.")
######################################################################################################################



#############################################Gene Ontology Enrichment Analysis######################################
print("Gene Ontology Enrichment Starts...")

#Enter the annotation db for the data set
##annotation_db <- readline(prompt = "Enter the annotation library ")
annotation_db <- "hgu133plus2.db"

#Enter the number of best p-values by user
##num <- readline(prompt = "Enter the best number of p-values for ontology ")
##numPvalue <- as.numeric(num)
numPvalue <- 10

moduleLabels_exper5 <- dynamicMods_exper5
#Finding the frequency of the genes per cluster
#It has two clumns:
#    #moduleLabels_STS: Numeric label of clusters
#    #Freq: Number of genes per cluster
clusters_label_freq_exper5 <- as.data.frame(table(moduleLabels_exper5))

#Creating the sumary data frame for best numPvalue ontology outcomes of the clusters
clusters_GO_exper5 <- data.frame(matrix(nrow = 0, ncol = 9))
summary_column_names <- c("clusterNum", "GOtype", "GOID" ,"Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term")
colnames(clusters_GO_exper5) <- summary_column_names

##########Calling gene ontology function per cluster
#Grouping the genes belong to the same cluster and then calling the gene ontology function
for(i in 1:dim(clusters_label_freq_exper5)[1]){
  #Extracting the i-th label
  label <- clusters_label_freq_exper5$moduleLabels_exper5[i]
  label <- droplevels(label)
  label <- as.numeric(levels(label))
  #
  ##label <- 11
  cluster <- -1
  for(j in 1:length(moduleLabels_exper5)){
    if(moduleLabels_exper5[j] == label){
      cluster = rbind(cluster, EntrezID_NoOutlier[j])
    }
  }
  cluster <- cluster[-1]
  if(length(cluster) != clusters_label_freq_exper5$Freq[i])
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
  clusters_GO_exper5 = rbind.fill(clusters_GO_exper5, res)
  
}

print("Gene Ontology Enrichment is done.")

######################################################################################################################




#############################################Saving The Results#######################################################
#Save the module colors and labels and the summary of annotation
caption_title <- paste("GSE", GSE, "_InnerProduct_Experiment5_3000.RData" ,sep = "")
save(moduleLabels_exper5, geneTree_exper5, clusters_GO_exper5, clusters_label_freq_exper5, file = caption_title )
#Followiing items will be saved:
#     moduleLabels_exper5           =>  numeric label of the clusters
#     geneTree_exper5               =>  denrogram of the clustering
#     clusters_GO_exper5            =>  Summary of the gene ontology annotation
#     clusters_label_freq_exper5    => table of the clusters label and their corresponding frequncies
######################################################################################################################
