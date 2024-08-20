###
geNorm <- function(df, verbose = T, PlotIt = TRUE)
{
  nStep <- length(colnames(df)) - 1
  
  geneOrder <- vector()
  AvgStability <- vector()
  
  for (iSTEP in 1:nStep)
  {
    matV <- calLogComp(df)
    matSDV <- calSDV(matV)
    
    eval <- matrix(unlist(strsplit(colnames(matV), split = "_", fixed = T)), ncol = 2, byrow = T)
    matM <- calM(df, dfeval = eval, mSDV = matSDV)
    
    if(verbose == TRUE){
      cat("##################################################################", "\n")
      cat("\n")
      cat("Step ", iSTEP, ":", "\n")
      cat("\n")
      cat("Stability Values:", "\n")
      print(matM[1,])
      cat("Average Stability:", mean(matM), "\n")
      cat("Lowest stability:", colnames(matM)[which.max(matM)], "\n")
      cat("\n")
    }
    
    geneOrder <- append(geneOrder, values = colnames(matM)[which.max(matM)])
    AvgStability <- append(AvgStability, values = mean(matM))
    
    df <- df[ , -which(colnames(df) %in% c(colnames(matM)[which.max(matM)]))]
    
  }
  
  geneOrder <- append(geneOrder, values = setdiff(colnames(matM), geneOrder))
  AvgStability <- append(AvgStability, values = tail(AvgStability, 1))
  dfStab <- data.frame(geneOrder, AvgStability)
  
  if (PlotIt == TRUE){
    plot(dfStab$AvgStability, pch = 21, bg = "gray50", ylab = "Average Expression Stability M", xaxt = "n", xlab = "")
    axis(side = 1, labels = dfStab$geneOrder, at = 1:length(dfStab$geneOrder), las = 2)
  }
  
  return(dfStab)
}

###
calLogComp <- function(df)
{
  matCT <- as.matrix(df)
  vV <- vector()
  vIJ <- vector()
  
  for (i in 1:ncol(matCT))
  {
    for (j in (i+1):ncol(matCT))
    {
      if (j > ncol(matCT)){break}
      for (iROW in 1:nrow(matCT))
      {
        V <- log2(matCT[iROW, i] / matCT[iROW, j])
        vV <- append(vV, values = V)
      }
      vIJ <- append(vIJ, values = paste(colnames(matCT)[i], colnames(matCT)[j], sep = "_"))
    }
  }
  
  matV <- matrix(data = vV, nrow = nrow(matCT), byrow = FALSE)
  colnames(matV) <- vIJ
  
  return(matV)
}

###
calSDV <- function(mV)
{
  vSDV <- vector()
  for (iCOL in 1:ncol(mV))
  {
    vSDV <- append(vSDV, values = sd(mV[,iCOL]))
  }
  
  matSDV <- matrix(data = vSDV, nrow = 1)
  colnames(matSDV) <- colnames(mV)
  return(matSDV)
}

###
calM <- function(df, dfeval, mSDV)
{
  vM <- vector()
  vGENE <- vector()
  for (iGENE in colnames(df))
  {
    iCOLS <- which(dfeval[,1] == iGENE)
    iCOLS <- append(iCOLS, values = which(dfeval[,2] == iGENE))
    
    vM <- append(vM, values = mean(mSDV[,iCOLS]))
    vGENE <- append(vGENE, values = iGENE)
  }
  
  matM <- matrix(data = vM, nrow = 1)
  colnames(matM) <- vGENE
  return(matM)
}