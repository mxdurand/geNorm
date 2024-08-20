# geNorm: Transcript accumulation normalization by housekeeping genes

Function `geNorm` takes data frame of standardized gene expression as input and returns gene stability.

Input data frame must have all genes in separate columns, and replicates as rows. Gene names as columns names, and replicate ID as row names. 

If verbose = TRUE, prints at each step information of stbility.
If PlotIt = TRUE, plot mean stability per gene. 
