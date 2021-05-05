##scalefree <- function(data, dataIsExpr <- TRUE, weights <- NULL, RsquaredCut <- 0.85,
##                      powerVector <- c(seq(1, 10, by = 1), seq(12, 20, by = 20)), removeFirst <- FALSE, nBreaks <- 10,
##                      blockSize <- NULL, corFnc <- cor, corOptions <- list(use = 'p'), networkType <- "unsigned",
##                      moreNetworkConcepts <- FALSE, gcInterval <- NULL, verbose <- 0, indent <- 0){
  #Here we find the power for which the adjacency matrix become scale free
  #
  #Arguments
  #data: Adjacency matrix of inner product (symmetric, squared, with diagonal equal to 1)
  #powerVector: integer vector of the powers 
  #
  #
  #Values
  #
  #
  #Niloofar Aghaieabiane
  #July 2019
  
  #initialization
  
  dataIsExpr <- TRUE
  weights <- NULL
  RsquaredCut <- 0.85
  #powerVector <- c(seq(1, 10, by = 1), seq(12, 20, by = 20))
  removeFirst <- FALSE
  nBreaks <- 10
  blockSize <- NULL
  corFnc <- cor
  corOptions <- list(use = 'p')
  networkType <- "unsigned"
  moreNetworkConcepts <- FALSE
  gcInterval <- NULL
  verbose <- 5
  indent <- 0

  #data <- matrix(c(1, 0.8, 0.3, 0.8,1,0.2,0.3,0.2,1), 3, 3, byrow = TRUE)
  #powerVector <- c(c(1:10))
 
  powerVector <- sort(powerVector)
  nGenes <- ncol(data)
  
  #check the adjacency matrix
  if(any(diag(data) != 1))
    stop("The diag of the adjacency matrix is not equal to 1")
  
  colname1 <- c("Power", "SFT.R.sq", "slope", "truncated R.sq", "mean(k)", "median(k)", "max(k)")
  
  datout <- data.frame(matrix(666, nrow = length(powerVector), ncol = length(colname1)))
  names(datout) <- colname1
  datout[, 1] <- powerVector
  
  datk <- matrix(0, nrow = nGenes, ncol = length(powerVector))
  
  nPowers <- length(powerVector)
  
  corx <- data
  corxPrev <- matrix(1, nrow = nrow(corx), ncol = ncol(corx))
  powerVector1 <- c(0, head(powerVector, -1))
  powerSteps <- powerVector - powerVector1
  uniquePowerSteps <- unique(powerSteps)
  corxPowers <- lapply(uniquePowerSteps, function(p) corx^p)
  names(corxPowers) <- uniquePowerSteps
  
  for(j in 1:nPowers){
    corxCur <- corxPrev * corxPowers[[as.character(powerSteps[j])]]
    datk[, j] <- colSums(corxCur, na.rm = TRUE) - 1
    corxPrev <- corxCur
  }
  
  for(i in c(1:length(powerVector))){
    khelp <- datk[, i]
    SFT1 <- scaleFreeFitIndex(k = khelp, nBreaks = nBreaks, removeFirst = removeFirst)
    datout[i,2] <- SFT1$Rsquared.SFT
    datout[i,3] <- SFT1$slope.SFT
    datout[i,4] <- SFT1$truncatedExponentialAdjRsquared
    datout[i,5] <- mean(khelp, na.rm = TRUE)
    datout[i,6] <- median(khelp, na.rm = TRUE)
    datout[i,7] <- max(khelp, na.rm = TRUE)
  }
  
  fitIndices <- data.frame(datout)
  sft <- list(fitIndices = data.frame(datout))
  ##print(signif(data.frame(datout), 3))
  ##ind1 = datout[, 2] > RsquaredCut
  ##indcut = NA
  ##indcut = if (sum(ind1) > 0) 
  ##  min(c(1:length(ind1))[ind1])
  ##        else indcut 
  ##powerEstimate <- powerVector[indcut][[1]]
  ##gc()
  ##list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
  
  
##}