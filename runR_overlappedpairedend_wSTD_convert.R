library(doParallel)

Community <- read.table("OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv", header=T, row.names=1)
Standard <- read.table("OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_standard.tsv", header=T, row.names=1)

stopifnot(nrow(Community) == nrow(Standard))

# Copy number per 1uL of internal standard
StandardCopy <- c(10, 20, 40, 80)

slope <- data.frame()
# Function definition
CalcSlope <- function(x){
  return(lm(as.numeric(Standard[x,]) ~ StandardCopy + 0)$coefficients)
}
ConvertReads <- function(x){
  conv <- Community[x,] / slope[x,1]
  return(conv)
}

# Leastsquare regression
slope <- data.frame()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
slope <- foreach(i = 1:nrow(Standard), .combine=rbind, .inorder=T, .export=c("Standard", "StandardCopy", "CalcSlope")) %dopar% {
	CalcSlope(i)
}
stopCluster(cl)

# Convert number of reads based on number of internal standard reads
Converted <- data.frame()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
Converted <- foreach(i = 1:nrow(Community), .combine=rbind, .inorder=T, .export=c("Community", "ConvertReads", "slope")) %dopar% {
	ConvertReads(i)
}
stopCluster(cl)

Converted[is.na(Converted)] <- NA
colnames(Converted) <- colnames(Community)
# Converted to DNA copy number per 1L (Filtered water amount = 1L, DNA extract = 200uL)
Converted <- Converted * 200 * (1 / 1)

# Save converted data
temp <- as.data.frame(row.names(Converted), row.names=row.names(Converted))
colnames(temp) <- "samplename"
write.table(cbind(temp, Converted), "OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_fishes_converted.tsv", sep="\t", append=F, quote=F, row.names=F, col.names=T, na="NA")
write.table(cbind(temp, slope), "OverlappedPairedEnd_wSTD_11_ClaidentResults/slope.tsv", sep="\t", append=F, quote=F, row.names=F, col.names=F, na="NA")
