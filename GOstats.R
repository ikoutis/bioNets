GO_analysis <- function(universe_geneIDs, cluster_geneIDs, annotation_db, nPVal, index, label, keepFiles = FALSE){

###universe_geneIDs <- EntrezID_NoOutlier
###cluster_geneIDs <- EntrezID_NoOutlier
###nPVal <- 10
###index <- 1
###label <- 1
###keepFiles <- FALSE
###annotation_db <- "hgu133plus2.db"

  #Here we perform the gene ontlogy for a given cluster
  #
  ##Arguments
  #   annotation_db: annotation library
  #   cluster_geneIDs: a set of gene EntrezIDs
  #   universe_geneIDs: all the gene EntrezIDs
  #   annotation_db: annotation library
  #   nPVal : number of best p-values to be returned
  #
  ##Values
  #
  #
  #
  #Niloofar Aghaieabiane
  #Jun 2019
  #
  #installing BiocManager
  
  #source("http://bioconductor.org/biocLite.R")
  #
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  #BiocManager::install()
  #
  if("GO.db" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("GO.db")}
  #biocLite("GO.db")
  library(GO.db)
  #
  if("annotate" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("annotate")}
  #biocLite("annotate")
  library(annotate)
  #
  if("genefilter" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("genefilter")}
  #biocLite("genefilter")
  library(genefilter)
  #
  if("GOstats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("GOstats")}
  #biocLite("GOstats")
  library(GOstats)
  #
  if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.package("RColorBrewer")}
  library(RColorBrewer)
  #
  if("xtable" %in% rownames(installed.packages()) == FALSE) {install.package("xtable")}
  library(xtable)
  #
  #if("Rgraphiviz" %in% rownames(installed.packages()) == FALSE) {install.package("Rgraphviz")}
  library(Rgraphviz)
  #
  #if("plyr" %in% rownames(installed.packages()) == FALSE) {install.package("plyr")}
  library(plyr)
  # 
  #if("BiocInstaller" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("BiocInstaller")}
  #biocLite("BiocInstaller")
  #library(BiocInstaller)
  #
  if(annotation_db %in% rownames(installed.packages()) == FALSE) {BiocManager::install(annotation_db)}
  #biocLite(annotation_db)
  library(annotation_db, character.only = TRUE)

  #biocLite("hgu133plus2.db")
  #BiocManager::install("hgu133plus2.db")
  #library("hgu133plus2.db")
  
  #Initialization
  test_direction <- c("over", "under")
  gene_ontology <- c("BP", "CC", "MF")

  #Creating the sumary data frame for all possible ontology outcomes 
  GO_summary <- data.frame(matrix(nrow = 0, ncol = 9))
  summary_column_names <- c("clusterNum", "GOtype", "GOID" ,"Pvalue", "OddsRatio", "ExpCount", "Count", "Size", "Term")
  colnames(GO_summary) <- summary_column_names
  
  
  #performing gene ontology for all six possible combination of diraction and gene ontology
  for(direction in test_direction){
    temp_direction <- direction
    #
    for(onto in gene_ontology){
      temp_ontology <- onto
      #
      hgCutoff <- 0.001
      params <- new("GOHyperGParams",
                    geneIds = cluster_geneIDs,
                    universeGeneIds = universe_geneIDs,
                    annotation = annotation_db,
                    ontology = temp_ontology,
                    #pvalueCutoff = hgCutoff,
                    conditional = FALSE,
                    testDirection = temp_direction)
      hg <- hyperGTest(params)
      if(keepFiles == TRUE)
        assign(paste(index, "_", temp_ontology,"_", temp_direction, sep = "" ), hg)
      #
      summary_hg <- summary(hg)
      if(dim(summary_hg)[1] != 0){
        summary_hg$clusterNum <- as.numeric(label)
        summary_hg$GOtype <- NaN
        summary_hg <- summary_hg[,c(8,9,1,2,3,4,5,6,7)]
        GO_type <- paste(direction, onto, sep = "")
        summary_hg$GOtype <- GO_type
        col_temp <- paste("GO", temp_ontology, "ID", sep = "")
        colnames(summary_hg)[colnames(summary_hg)== col_temp] <- "GOID"
        #
        GO_summary <- rbind.fill(GO_summary,summary_hg)
        rm(list = c("summary_hg")) 
      }
    }
  }
  hg
  #sort the gene ontolgy result of all combination based on the p-values
  GO_summary <- GO_summary[order(GO_summary$Pvalue),]
  
  #pick the numPVal given by user
  GO_result <- GO_summary[1:nPVal, ]
  newlist <- list("hg" = hg, "GO_result" = GO_result)
  return(GO_result)
  #return(newlist)
}


