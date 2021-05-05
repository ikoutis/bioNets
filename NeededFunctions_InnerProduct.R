sigOnObservation <- function(dataset){
  ##Taking the sigmoid function on the observation 
  #
  #Argument
  # An expression matrix of #samples by #genes
  #
  #Values
  # A matrix of the same dimension of the input that sigmoid function is applied on it
  #
  # Niloofar Aghaieabiane
  # July 2019
  
  #loading required libraries and packages
  library(pracma)
  
  print("Calculating the sigmoid function for each entry in expression data set...")
  sig_dataset <- matrix(NaN, dim(dataset)[1], dim(dataset)[2])
  #
  for(j in 1:(dim(dataset)[2])){
    #claculating variance per gene (i.e. column)
    variance <- var(dataset[,j])
    #calculating mean per gene (i.e. column)
    average <- mean(dataset[,j])
    #
    for(i in 1:dim(dataset)[1]){
      sig_dataset[i,j] <- sigmoid(dataset[i,j], a = (1/variance), b = average )
    }
  }
  
  if(dim(dataset)[1]==dim(sig_dataset)[1] & dim(dataset)[2]==dim(sig_dataset)[2])
    print("The calculation of sigmoid function is done.")
  
  return(sig_dataset)
  
}

standardization <- function(dataset){
  ## Here we apply the z-tranform the enetry of the dataset (i.e. (x - mean)/sd )
  #The dataset must be symmetric and squared before and after z-tranfrom
  #
  #Argument
  # A symmetric squared #gene by #gene matrix
  #
  #Values
  #z-transform of the matrix of #gene by #gene (must be symmetric and squared)
  
  stand_matrix <- matrix(NaN, dim(dataset)[1], dim(dataset)[2])
  
  for(i in 1:dim(dataset)[1]){
    #
    average <- mean(dataset[i,])
    standard_dev <- sd(dataset[i,])
    for(j in i:dim(dataset)[2]){
      if(i != j){
        stand_matrix[i,j] <- (dataset[i,j] - average)/standard_dev
        stand_matrix[j,i] <- (dataset[i,j] - average)/standard_dev 
      }
    }
  }
  diag(stand_matrix) <- 1
  
  if(sum(is.na(stand_matrix)) != 0)
    stop("error, there is a NaN value in stand_matrix")
  
  if(isSymmetric(stand_matrix) == TRUE)
    print("Standardization is done.")
  
  return(stand_matrix)
}

normalization <- function(dataset){
  ## Here we normalize the enetry of the dataset (i.e. (x - x_min/x_max - x_min )
  #The dataset must be symmetric and squared before and after normalizing
  #
  #Argument
  # A symmetric squared #gene by #gene matrix
  #
  #Values
  #Normalized matrix of #gene by #gene (must be symmetric and squared)
  
  norm_matrix <- matrix(NaN, dim(dataset)[1], dim(dataset)[2])
  
  for(i in 1:dim(dataset)[1]){
    #
    minimum <- min(dataset[i,])
    maximum <- max(dataset[i,])
    for(j in i:dim(dataset)[2]){
      if(i != j){
        norm_matrix[i,j] <- (dataset[i,j] - minimum)/(maximum - minimum)
        norm_matrix[j,i] <- (dataset[i,j] - minimum)/(maximum - minimum) 
      }
    }
  }
  diag(norm_matrix) <- 1
  
  if(sum(is.na(norm_matrix)) != 0)
    stop("error, there is a NaN value in norm_matrix")
  
  if(isSymmetric(norm_matrix) == TRUE)
    print("Normalization is done.")
  
  return(norm_matrix)
}


maxDivide <- function(dataset){
  ## Here we divide the enetry of the dataset by the maximum element (i.e. (x - x_min/(x_max - x_min) )
  #The dataset must be symmetric and squared before and after normalizing
  #
  #Argument
  # A symmetric squared #gene by #gene matrix
  #
  #Values
  #A divide by maximum matrix of #gene by #gene (must be symmetric and squared) 
  #so that all the elements are between 0 and 1
  
  #maxDivide_matrix <- matrix(NaN, dim(dataset)[1], dim(dataset)[2])
  
  
  maximum <- max(dataset)
  maxDivide_matrix <- dataset /max(dataset)
  
  diag(maxDivide_matrix) <- 1
  
  if(sum(is.na(maxDivide_matrix)) != 0)
    stop("error, there is a NaN value in maxDivide_matrix")
  
  if(isSymmetric(maxDivide_matrix) == TRUE)
    print("Divide by maximum is done.")
  
  return(maxDivide_matrix)
}