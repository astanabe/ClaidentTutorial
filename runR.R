library(tidyverse)
library(ggsci)
library(foreach)
library(doParallel)
library(vegan)
library(mpmcorrelogram)
library(geosphere)
library(scales)

# Make output directory
dir.create("11_RAnalysisResults")

# Make species-level Barplot
pdf("11_RAnalysisResults/top50species.pdf", width=14, height=10)
top50species <- read.table("10_ClaidentResults/sample_top50species_matrix_fishes.tsv", header=T)
temp <- ggplot(top50species, aes(x=samplename, y=nreads, fill=fct_rev(species)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="species")
temp <- temp + theme_test()
temp <- temp + theme(axis.text=element_text(angle=90))
plot(temp)
dev.off()

# Make family-level Barplot
pdf("11_RAnalysisResults/top50family.pdf", width=14, height=10)
top50family <- read.table("10_ClaidentResults/sample_top50family_matrix_fishes.tsv", header=T)
temp <- ggplot(top50family, aes(x=samplename, y=nreads, fill=fct_rev(family)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="family")
temp <- temp + theme_test()
temp <- temp + theme(axis.text=element_text(angle=90))
plot(temp)
dev.off()

# Read community data matrix
Community <- read.table("10_ClaidentResults/sample_species_matrix_fishes.tsv", header=T, row.names=1)

# Read metadata
Metadata <- read.table("Metadata.tsv", header=T, row.names=1)

# Draw species accumulation curve
SpecAccum <- specaccum(Community)
pdf("11_RAnalysisResults/specaccum.pdf", width=7, height=7)
plot(SpecAccum, xlab="number of samples", ylab="number of species", main="species accumulation curve")
dev.off()

# Draw rarefaction curves
pdf("11_RAnalysisResults/rarecurve.pdf")
rarecurve(Community, step=10, xlab="number of seqs", ylab="number of species", main="rarefaction curves")
dev.off()

# Coverage-based rarefaction
## make rareslopelist using all cpu cores
rareslopelist <- list()
cl <- makeCluster(detectCores())
registerDoParallel(cl)
rareslopelist <- foreach(i = 1:nrow(Community), .packages="vegan") %dopar% {
	rareslope(Community[i,], seq(1, (sum(Community[i,]) - 1), by=1))
}
stopCluster(cl)
## find minimum coverage sample
getmincov <- c()
for(i in 1:nrow(Community)) {
    getmincov[i] <- rareslopelist[[i]][length(rareslopelist[[i]])]
}
## echo minimum coverage
(1 - max(getmincov)) * 100
## set target slope
## to demonstrate coverage-based rarefaction, cvr is set to 0.05 (95% coverage)
cvr <- 0.05
## in the actual analysis, the following value is recommended
#cvr <- max(getmincov)
## define function
cvrfun <- function(x) {min(which(x <= cvr))}
## get number of seqs of target coverage
cvrrare <- unlist(lapply(rareslopelist, cvrfun))
# make rarefied community data
RarefiedCommunity <- rrarefy(Community, cvrrare)
write.table(RarefiedCommunity, "11_RAnalysisResults/RarefiedCommunity.tsv", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

# Make binary community data
BinaryRarefiedCommunity <- data.frame()
BinaryRarefiedCommunity <- replace(RarefiedCommunity, RarefiedCommunity > 0, 1)
write.table(BinaryRarefiedCommunity, "11_RAnalysisResults/BinaryRarefiedCommunity.tsv", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

# Make Bray-Curtis distance matrix
BrayCurtis <- vegdist(RarefiedCommunity, method="bray")

# Make Jaccard distance matrix
Jaccard <- vegdist(RarefiedCommunity, method="jaccard")

# Make Jaccard distance matrix using binary data
BinaryJaccard <- vegdist(BinaryRarefiedCommunity, method="jaccard")

# Make Raup-Crick distance matrix using binary data
BinaryRaupCrick <- raupcrick(BinaryRarefiedCommunity, null="r1", nsimul=999)

# PERMANOVA
sink("11_RAnalysisResults/BrayCurtisPERMANOVA.txt", split=T)
adonis(BrayCurtis ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("11_RAnalysisResults/JaccardPERMANOVA.txt", split=T)
adonis(Jaccard ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("11_RAnalysisResults/BinaryJaccardPERMANOVA.txt", split=T)
adonis(BinaryJaccard ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("11_RAnalysisResults/BinaryRaupCrickPERMANOVA.txt", split=T)
adonis(BinaryRaupCrick ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()

# NMDS
## fit NMDS
BrayCurtisNMDS <- metaMDS(BrayCurtis, k=2, trymax=100)
JaccardNMDS <- metaMDS(Jaccard, k=2, trymax=100)
BinaryJaccardNMDS <- metaMDS(BinaryJaccard, k=2, trymax=100)
BinaryRaupCrickNMDS <- metaMDS(BinaryRaupCrick, k=2, trymax=100)
## renew NMDS by better-fitted results
for(i in 1:100) {
  set.seed(i)
  temp <- metaMDS(BrayCurtis, k=2, trymax=100)
  if(temp$stress < BrayCurtisNMDS$stress) {
    BrayCurtisNMDS <- temp
  }
}
for(i in 1:100) {
  set.seed(i)
  temp <- metaMDS(Jaccard, k=2, trymax=100)
  if(temp$stress < JaccardNMDS$stress) {
    JaccardNMDS <- temp
  }
}
for(i in 1:100) {
  set.seed(i)
  temp <- metaMDS(BinaryJaccard, k=2, trymax=100)
  if(temp$stress < BinaryJaccardNMDS$stress) {
    BinaryJaccardNMDS <- temp
  }
}
for(i in 1:100) {
  set.seed(i)
  temp <- metaMDS(BinaryRaupCrick, k=2, trymax=100)
  if(temp$stress < BinaryRaupCrickNMDS$stress) {
    BinaryRaupCrickNMDS <- temp
  }
}
## fit environmental data
BrayCurtisNMDSenv <- envfit(BrayCurtisNMDS, Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
JaccardNMDSenv <- envfit(JaccardNMDS, Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
BinaryJaccardNMDSenv <- envfit(BinaryJaccardNMDS, Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
BinaryRaupCrickNMDSenv <- envfit(BinaryRaupCrickNMDS, Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
## draw NMDS
pdf("11_RAnalysisResults/NMDS.pdf")
ordiplot(BrayCurtisNMDS, type="n")
orditorp(BrayCurtisNMDS, display="sites", air=0.1, cex=1)
plot(BrayCurtisNMDSenv, p.max=0.05)
ordiplot(JaccardNMDS, type="n")
orditorp(JaccardNMDS, display="sites", air=0.1, cex=1)
plot(JaccardNMDSenv, p.max=0.05)
ordiplot(BinaryJaccardNMDS, type="n")
orditorp(BinaryJaccardNMDS, display="sites", air=0.1, cex=1)
plot(BinaryJaccardNMDSenv, p.max=0.05)
ordiplot(BinaryRaupCrickNMDS, type="n")
orditorp(BinaryRaupCrickNMDS, display="sites", air=0.1, cex=1)
plot(BinaryRaupCrickNMDSenv, p.max=0.05)
dev.off()

# Spatial Mantel correlogram analysis
## make geographical distance matrix
Geodist <- distm(cbind(Metadata$Longitude, Metadata$Latitude), fun=distGeo)
## run analysis
BrayCurtisGeoMCA <- mpmcorrelogram(BrayCurtis, Geodist, method="spearman", permutations=999)
JaccardGeoMCA <- mpmcorrelogram(Jaccard, Geodist, method="spearman", permutations=999)
BinaryJaccardGeoMCA <- mpmcorrelogram(BinaryJaccard, Geodist, method="spearman", permutations=999)
BinaryRaupCrickGeoMCA <- mpmcorrelogram(BinaryRaupCrick, Geodist, method="spearman", permutations=999)
## draw analysis results
xval <- c()
pdf("11_RAnalysisResults/GeoMCA.pdf")
tipos <- BrayCurtisGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BrayCurtisGeoMCA$breaks) - 1)) {
    xval[i] <- (BrayCurtisGeoMCA$breaks[i] + BrayCurtisGeoMCA$breaks[i + 1]) / 2
}
plot(xval, BrayCurtisGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisGeoMCA$breaks, n=1)), xlab="geographic distance", ylab="Mantel correlation", main="BrayCurtisGeoMCA")
abline(v=BrayCurtisGeoMCA$breaks)
tipos <- JaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(JaccardGeoMCA$breaks) - 1)) {
    xval[i] <- (JaccardGeoMCA$breaks[i] + JaccardGeoMCA$breaks[i + 1]) / 2
}
plot(xval, JaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardGeoMCA$breaks, n=1)), xlab="geographic distance", ylab="Mantel correlation", main="JaccardGeoMCA")
abline(v=JaccardGeoMCA$breaks)
tipos <- BinaryJaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryJaccardGeoMCA$breaks) - 1)) {
    xval[i] <- (BinaryJaccardGeoMCA$breaks[i] + BinaryJaccardGeoMCA$breaks[i + 1]) / 2
}
plot(xval, BinaryJaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardGeoMCA$breaks, n=1)), xlab="geographic distance", ylab="Mantel correlation", main="BinaryJaccardGeoMCA")
abline(v=BinaryJaccardGeoMCA$breaks)
tipos <- BinaryRaupCrickGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryRaupCrickGeoMCA$breaks) - 1)) {
    xval[i] <- (BinaryRaupCrickGeoMCA$breaks[i] + BinaryRaupCrickGeoMCA$breaks[i + 1]) / 2
}
plot(xval, BinaryRaupCrickGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickGeoMCA$breaks, n=1)), xlab="geographic distance", ylab="Mantel correlation", main="BinaryRaupCrickGeoMCA")
abline(v=BinaryRaupCrickGeoMCA$breaks)
dev.off()

# Temporal Mantel correlogram analysis
## make date distance matrix
Datedist <- dist(as.Date(Metadata$Date), diag=T, upper=T)
## run analysis
BrayCurtisDateMCA <- mpmcorrelogram(BrayCurtis, Datedist, method="spearman", permutations=999)
JaccardDateMCA <- mpmcorrelogram(Jaccard, Datedist, method="spearman", permutations=999)
BinaryJaccardDateMCA <- mpmcorrelogram(BinaryJaccard, Datedist, method="spearman", permutations=999)
BinaryRaupCrickDateMCA <- mpmcorrelogram(BinaryRaupCrick, Datedist, method="spearman", permutations=999)
## draw analysis results
xval <- c()
pdf("11_RAnalysisResults/DateMCA.pdf")
tipos <- BrayCurtisDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BrayCurtisDateMCA$breaks) - 1)) {
    xval[i] <- (BrayCurtisDateMCA$breaks[i] + BrayCurtisDateMCA$breaks[i + 1]) / 2
}
plot(xval, BrayCurtisDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisDateMCA$breaks, n=1)), xlab="date interval", ylab="Mantel correlation", main="BrayCurtisDateMCA")
abline(v=BrayCurtisDateMCA$breaks)
tipos <- JaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(JaccardDateMCA$breaks) - 1)) {
    xval[i] <- (JaccardDateMCA$breaks[i] + JaccardDateMCA$breaks[i + 1]) / 2
}
plot(xval, JaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardDateMCA$breaks, n=1)), xlab="date interval", ylab="Mantel correlation", main="JaccardDateMCA")
abline(v=JaccardDateMCA$breaks)
tipos <- BinaryJaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryJaccardDateMCA$breaks) - 1)) {
    xval[i] <- (BinaryJaccardDateMCA$breaks[i] + BinaryJaccardDateMCA$breaks[i + 1]) / 2
}
plot(xval, BinaryJaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardDateMCA$breaks, n=1)), xlab="date interval", ylab="Mantel correlation", main="BinaryJaccardDateMCA")
abline(v=BinaryJaccardDateMCA$breaks)
tipos <- BinaryRaupCrickDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryRaupCrickDateMCA$breaks) - 1)) {
    xval[i] <- (BinaryRaupCrickDateMCA$breaks[i] + BinaryRaupCrickDateMCA$breaks[i + 1]) / 2
}
plot(xval, BinaryRaupCrickDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickDateMCA$breaks, n=1)), xlab="date interval", ylab="Mantel correlation", main="BinaryRaupCrickDateMCA")
abline(v=BinaryRaupCrickDateMCA$breaks)
dev.off()
