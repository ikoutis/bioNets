# A novel calibration step in gene co-expression network construction
A novel step toward the gene co-expression network (GCN) construction. 

## Introduction
A new step toward the gene co-expression network is introduced, that it makes the use of Sigmoid function of raw data of gene expression.
This step in effect calibrates the expression levels of individual genes, before computing pairwise similarities. Further, the similarity is computed as an inner-product of positive vectors.
In multiple experiments, this seems to provide a significant improvement, as measured by aggregate p-values of the gene ontology term enrichment of the computed modules.
The prevously well-known framework for constructing GCNs, WGCNA, has been compared with three new framework based on this preprocessing step; Alpha, Beta, and Gamma.

The files corresponding with each framework is follow:  
Experiment_WGCNA.R, Experiment_Alpha.R, Experiment_Beta.R, Experiment_Gamma.R

### Experiment_WGCNA.R includes the following steps:
i. Calculate Pearson correlation on gene expression.
ii. Convert negative values into positive.
iii. Power the network (element wise) so that the underlying network becomes scale-free.
iv. Add topological overlap matrix (TOM) to the network.
v. Perform clustering.
vi. Perform GO. 

### Experiment_Alpha.R includes the following steps:
i. Apply the novel preprocessing step on raw data of gene expression. 
ii. Compute the inner product.
iii. Normalize the matrix from previouse step.
iv. Power the network so that it becomes scale-free.
v. Perform clustering.
vi. Perform GO. 

### Experiment_Beta.R includes the following steps:
i. Apply the novel preprocessing step on raw data of gene expression.
ii. Divide each vector of gene expression by its norm.
iii. Compute the inner product.  
iv. Power the network so that it becomes scale-free.
v. Perform clustering.
vi. Perform GO. 

### Experiment_Gamma.R includes the following stes:
i. Follow the steps of Experiment_Beta.R.
ii. Add TOM to the network.

For all the four frameworks the following setting is the same:
* function "cutreeDynamics" with setting method = "tree", and minModuleSize = 30 from package "dynamicTreeCut.R" has been used for clustering.
* function "hyperGTest" from package "GOstats.R" with all the posibilities of "underBP", "overBP", "underCC", "overCC", "underMF", "overMF" has been used for gene ontology.


### Input Data (for all the frameworks)
A CSV file in which each entry i and j indicates the expression gene j in sample i.
A CSV file indicates the gene EntrezID corresponding to the prevouse CSV file.

### Output Data (for all the frameworks)
A data frame (*.RData file) containing the gene ontology (GO) information of the clusters found by the framework;
ClusterNum, GOtype, GOID, Pvalue, OddsRatio, ExpCount, Count, Size, and Term.

### Link to the expression data sets:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34400
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129166
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27948
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137394
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30140
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27211


 @na396


