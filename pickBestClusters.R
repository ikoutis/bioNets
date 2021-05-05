pickBestClusters_2 <- function(df, index, numClusters, Method, x_axis_name){
  #This function firstly calculate the geometric mean for the index number of p-values per cluaster
  # and then it returns the numClsuters geo mean of pvalues among the all the clusters that are the best
  #
  #For example, if method produces 20 clusters and each cluster has 12 P-values, and the index is 4 and the numClusters is 17, then
  # it firstly pick the 4 best p-values for each cluster
  # calculate the geometric mean over the 4 best Pvalues for each cluster
  # among the geometric mean for 20 clusters, it return the 17 best geometric mean 
  #
  #Argument
  # df: dataframe of GO results of the different experiments
  # index: number of best p-value that should be considered for each cluster
  # numClusters: nuumber of best clusters that has been retuned
  #
  #Values
  #out: a dataframe contaning the clusterNum, and Pvalue (geometric mean over the index number of best p-values for each cluster)
  #
  
  
  
  
  ###df <- clusters_GO_TOM
  num <- index
  
  clusters_label <- unique(df$clusterNum)
  
  #defining the output for the df
  temp_out <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(temp_out) <- c("clusterNum", "Pvalue")
  
  for(i in 1:length(clusters_label) ){
    label <- clusters_label[i]
    if(label != 0){
      temp <- df[c(df$clusterNum == label),]
      temp <- temp[order(temp$Pvalue),]
      temp <- temp[1:num,]
      keep_columns <- c("clusterNum", "Pvalue")
      temp <- temp[keep_columns]
      temp$Pvalue <- log10(t(t(temp$Pvalue)))
      temp$Pvalue <- -1 * temp$Pvalue
      
      temp_out <- rbind(temp_out, temp)
      rm(temp)
    }
  }
  out <- aggregate(Pvalue ~ clusterNum, temp_out, mean)
  colnames(out) <- c("clusterNum","geo_mean_Pvalue")
  
  ##Sorting the out_TOM data frame in respect of the Pvalue column to be able to pick the 10 clusters with highest quality
  out <- out[order(out$geo_mean_Pvalue, decreasing = TRUE),]
  out <- out[1:numClusters,]
  out$Method <- Method
  out$x_axis <- x_axis_name 
  
  return(out)
}




NumOfGenesInCluster <- function(df1, df2, column_name, numClusters){
  #This function cacultes the number of genes for each cluster
  #
  ##Argument:
  # df1: a dataframe retuned by pickBestClusters function(has four column: clusterNum, geo_mean_Pvalue, Method, x_axis)
  # df2: a dataframe of containing the clusters' label and the number of associated genes (has two columns: moduleLabels_method, and Freq)
  # column_name: name of the first column in the df2, either moduleLabels_TOM , moduleLabels_exper5 , moduleLabels_exper6
  #            , moduleLabels_exper11, moduleLabels_exper12
  #
  ##Values
  #
  
  df1$clusterSize <- NA
  
  for(i in 1:numClusters){
    
    cluster_name <- df1$clusterNum[i]
    temp <- df2[df2[[column_name]] == cluster_name, ]
    df1$clusterSize[i] <- temp$Freq
    
  }
}


