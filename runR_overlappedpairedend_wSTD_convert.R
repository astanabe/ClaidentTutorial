Community <- read.table("OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_fishes.tsv", header=T, row.names=1)
Standard <- read.table("OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_standard.tsv", header=T, row.names=1)

# Copy number per 1uL of internal standard
StandardCopy <- c(10, 20, 40, 80)

# Function definition
ConvertReads <- function(x){
  slope <- lm(as.numeric(Standard[x,]) ~ StandardCopy + 0)$coefficients
  conv <- Community[x,] / slope
  return(conv)
}

# Convert number of reads based on number of internal standard reads
Converted <- data.frame()
for(i in 1:nrow(Community)){
  Converted <- rbind(Converted, ConvertReads(i))
}
Converted[is.na(Converted)] <- NA
colnames(Converted) <- colnames(Community)
# Converted to DNA copy number per 1L (Filtered water amount = 1L, DNA extract = 200uL)
Converted <- Converted * 200 * (1 / 1)

# Save converted data
temp <- as.data.frame(row.names(Converted), row.names=row.names(Converted))
colnames(temp) <- "samplename"
write.table(cbind(temp, Converted), "OverlappedPairedEnd_wSTD_11_ClaidentResults/sample_otu_matrix_fishes_converted.tsv", sep="\t", append=F, quote=F, row.names=F, col.names=T, na="NA")
