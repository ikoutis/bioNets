Pick_Pvalues <- function(df, num){
# This function gives a data frame of gene ontology of the clusters and 
# it returns the num best p-values for each cluster and   
# it returns the sorted data frame of in terms of the p-values of per clusters 
#
# Arguments
# df: a dataframe of the gene ontology of the clusters
# num: number of the best pvalues per cluster
#
# Values
# best_pvalues: the num best p-values per cluster  
# df: A dataframe in which clusters are sorted based on their P-values



clusters_label <- unique(df$clusterNum)


best_pvalues <- -1

for(label in clusters_label){
  temp <- df[c(df$clusterNum == label),]
  temp <- temp[order(temp$Pvalue),]
  if( label != 0)
    best_pvalues <- rbind(best_pvalues, temp$Pvalue[1:num])
  df[c(df$clusterNum == label),] <- temp
}

df

best_pvalues <- best_pvalues[-1]
return(best_pvalues)

}
