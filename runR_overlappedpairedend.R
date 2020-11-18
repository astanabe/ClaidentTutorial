library(tidyverse)
library(ggsci)
library(foreach)
library(doParallel)
library(vegan)
library(mpmcorrelogram)
library(geosphere)
library(scales)

# Make output directory
dir.create("OverlappedPairedEnd_11_RAnalysisResults")

# Make species-level barplot
pdf("OverlappedPairedEnd_11_RAnalysisResults/barplottop50species.pdf", width=14, height=10)
top50species <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_top50species_nreads_fishes.tsv", header=T)
temp <- ggplot(top50species, aes(x=samplename, y=nreads, fill=fct_rev(species)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="species")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
plot(temp)
dev.off()

# Make family-level barplot
pdf("OverlappedPairedEnd_11_RAnalysisResults/barplottop50family.pdf", width=13, height=10)
top50family <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_top50family_nreads_fishes.tsv", header=T)
temp <- ggplot(top50family, aes(x=samplename, y=nreads, fill=fct_rev(family)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="family")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
plot(temp)
dev.off()

# Make species-level heatmap
pdf("OverlappedPairedEnd_11_RAnalysisResults/heatmapspecies.pdf", width=22, height=10)
commspecies <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_species_nreads_fishes.tsv", header=T)
commspecies$nreads[(commspecies$nreads == 0)] <- NA
temp <- ggplot(commspecies, aes(x=species, y=samplename, fill=nreads))
temp <- temp + geom_tile()
temp <- temp + scale_fill_gsea(na.value="white")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
plot(temp)
dev.off()

# Make family-level heatmap
pdf("OverlappedPairedEnd_11_RAnalysisResults/heatmapfamily.pdf", width=16, height=10)
commfamily <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_family_nreads_fishes.tsv", header=T)
commfamily$nreads[(commfamily$nreads == 0)] <- NA
temp <- ggplot(commfamily, aes(x=family, y=samplename, fill=nreads))
temp <- temp + geom_tile()
temp <- temp + scale_fill_gsea(na.value="white")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
plot(temp)
dev.off()

# Read community data matrix
Community <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_species_matrix_fishes.tsv", header=T, row.names=1)

# Draw species accumulation curve
SpecAccum <- specaccum(Community)
pdf("OverlappedPairedEnd_11_RAnalysisResults/specaccum.pdf", width=7, height=7)
plot(SpecAccum, xlab="number of samples", ylab="number of species", main="species accumulation curve")
dev.off()

# Draw rarefaction curves
pdf("OverlappedPairedEnd_11_RAnalysisResults/rarecurve.pdf", width=7, height=7)
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
## To demonstrate coverage-based rarefaction, cvr is set to 0.05 (95% coverage), but this is inappropreate in this case.
cvr <- 0.05
## In the actual analysis, the following value is recommended.
#cvr <- max(getmincov)
## define function
cvrfun <- function(x) {min(which(x <= cvr)) + 1}
## get number of seqs of target coverage
cvrrare <- unlist(lapply(rareslopelist, cvrfun))
# make rarefied community data
RarefiedCommunity <- rrarefy(Community, cvrrare)
write.table(RarefiedCommunity, "OverlappedPairedEnd_11_RAnalysisResults/RarefiedCommunity.tsv", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

# Make binary community data
BinaryRarefiedCommunity <- data.frame()
BinaryRarefiedCommunity <- replace(RarefiedCommunity, RarefiedCommunity > 0, 1)
write.table(BinaryRarefiedCommunity, "OverlappedPairedEnd_11_RAnalysisResults/BinaryRarefiedCommunity.tsv", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

# Make Bray-Curtis distance matrix
BrayCurtis <- vegdist(RarefiedCommunity, method="bray")

# Make Jaccard distance matrix
Jaccard <- vegdist(RarefiedCommunity, method="jaccard")

# Make Jaccard distance matrix using binary data
BinaryJaccard <- vegdist(BinaryRarefiedCommunity, method="jaccard")

# Make Raup-Crick distance matrix using binary data
BinaryRaupCrick <- as.dist(raupcrick(BinaryRarefiedCommunity, null="r1", nsimul=999))

# Read metadata
Metadata <- read.table("Metadata.tsv", header=T, row.names=1)

# PERMANOVA
sink("OverlappedPairedEnd_11_RAnalysisResults/BrayCurtisPERMANOVA.txt", split=T)
adonis(BrayCurtis ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("OverlappedPairedEnd_11_RAnalysisResults/JaccardPERMANOVA.txt", split=T)
adonis(Jaccard ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("OverlappedPairedEnd_11_RAnalysisResults/BinaryJaccardPERMANOVA.txt", split=T)
adonis(BinaryJaccard ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores())
sink()
sink("OverlappedPairedEnd_11_RAnalysisResults/BinaryRaupCrickPERMANOVA.txt", split=T)
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
pdf("OverlappedPairedEnd_11_RAnalysisResults/NMDS.pdf", width=7, height=7)
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
Geodist <- as.dist(distm(cbind(Metadata$Longitude, Metadata$Latitude), fun=distGeo))
## run analysis
BrayCurtisGeoMCA <- mpmcorrelogram(BrayCurtis, Geodist, method="spearman", permutations=999)
JaccardGeoMCA <- mpmcorrelogram(Jaccard, Geodist, method="spearman", permutations=999)
BinaryJaccardGeoMCA <- mpmcorrelogram(BinaryJaccard, Geodist, method="spearman", permutations=999)
BinaryRaupCrickGeoMCA <- mpmcorrelogram(BinaryRaupCrick, Geodist, method="spearman", permutations=999)
## draw analysis results
xval <- c()
pdf("OverlappedPairedEnd_11_RAnalysisResults/GeoMCA.pdf", width=7, height=7)
tipos <- BrayCurtisGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BrayCurtisGeoMCA$breaks) - 1)) {
    xval[i] <- (BrayCurtisGeoMCA$breaks[i] + BrayCurtisGeoMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Geodist, BrayCurtis, pch=16, cex=1, type="p", xlim=c(0, tail(BrayCurtisGeoMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Bray-Curtis distance", side=4, line=2.8)
par(new=T)
plot(xval, BrayCurtisGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BrayCurtisGeoMCA")
abline(v=BrayCurtisGeoMCA$breaks, lty=2)
tipos <- JaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(JaccardGeoMCA$breaks) - 1)) {
    xval[i] <- (JaccardGeoMCA$breaks[i] + JaccardGeoMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Geodist, Jaccard, pch=16, cex=1, type="p", xlim=c(0, tail(JaccardGeoMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Jaccard distance", side=4, line=2.8)
par(new=T)
plot(xval, JaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="JaccardGeoMCA")
abline(v=JaccardGeoMCA$breaks, lty=2)
tipos <- BinaryJaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryJaccardGeoMCA$breaks) - 1)) {
    xval[i] <- (BinaryJaccardGeoMCA$breaks[i] + BinaryJaccardGeoMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Geodist, BinaryJaccard, pch=16, cex=1, type="p", xlim=c(0, tail(BinaryJaccardGeoMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Jaccard distance of binary-transformed data", side=4, line=2.8)
par(new=T)
plot(xval, BinaryJaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BinaryJaccardGeoMCA")
abline(v=BinaryJaccardGeoMCA$breaks, lty=2)
tipos <- BinaryRaupCrickGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryRaupCrickGeoMCA$breaks) - 1)) {
    xval[i] <- (BinaryRaupCrickGeoMCA$breaks[i] + BinaryRaupCrickGeoMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Geodist, BinaryRaupCrick, pch=16, cex=1, type="p", xlim=c(0, tail(BinaryRaupCrickGeoMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Raup-Crick distance", side=4, line=2.8)
par(new=T)
plot(xval, BinaryRaupCrickGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BinaryRaupCrickGeoMCA")
abline(v=BinaryRaupCrickGeoMCA$breaks, lty=2)
dev.off()

# Temporal Mantel correlogram analysis
## make date distance matrix
Datedist <- dist(as.Date(Metadata$Date))
## run analysis
BrayCurtisDateMCA <- mpmcorrelogram(BrayCurtis, Datedist, method="spearman", permutations=999)
JaccardDateMCA <- mpmcorrelogram(Jaccard, Datedist, method="spearman", permutations=999)
BinaryJaccardDateMCA <- mpmcorrelogram(BinaryJaccard, Datedist, method="spearman", permutations=999)
BinaryRaupCrickDateMCA <- mpmcorrelogram(BinaryRaupCrick, Datedist, method="spearman", permutations=999)
## draw analysis results
xval <- c()
pdf("OverlappedPairedEnd_11_RAnalysisResults/DateMCA.pdf", width=7, height=7)
tipos <- BrayCurtisDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BrayCurtisDateMCA$breaks) - 1)) {
    xval[i] <- (BrayCurtisDateMCA$breaks[i] + BrayCurtisDateMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Datedist, BrayCurtis, pch=16, cex=1, type="p", xlim=c(0, tail(BrayCurtisDateMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Bray-Curtis distance", side=4, line=2.8)
par(new=T)
plot(xval, BrayCurtisDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisDateMCA$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BrayCurtisDateMCA")
abline(v=BrayCurtisDateMCA$breaks, lty=2)
tipos <- JaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(JaccardDateMCA$breaks) - 1)) {
    xval[i] <- (JaccardDateMCA$breaks[i] + JaccardDateMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Datedist, Jaccard, pch=16, cex=1, type="p", xlim=c(0, tail(JaccardDateMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Jaccard distance", side=4, line=2.8)
par(new=T)
plot(xval, JaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardDateMCA$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="JaccardDateMCA")
abline(v=JaccardDateMCA$breaks, lty=2)
tipos <- BinaryJaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryJaccardDateMCA$breaks) - 1)) {
    xval[i] <- (BinaryJaccardDateMCA$breaks[i] + BinaryJaccardDateMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Datedist, BinaryJaccard, pch=16, cex=1, type="p", xlim=c(0, tail(BinaryJaccardDateMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Jaccard distance of binary-transformed data", side=4, line=2.8)
par(new=T)
plot(xval, BinaryJaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardDateMCA$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BinaryJaccardDateMCA")
abline(v=BinaryJaccardDateMCA$breaks, lty=2)
tipos <- BinaryRaupCrickDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
for(i in 1:(length(BinaryRaupCrickDateMCA$breaks) - 1)) {
    xval[i] <- (BinaryRaupCrickDateMCA$breaks[i] + BinaryRaupCrickDateMCA$breaks[i + 1]) / 2
}
par(mar=c(5, 4, 4, 4))
plot(Datedist, BinaryRaupCrick, pch=16, cex=1, type="p", xlim=c(0, tail(BinaryRaupCrickDateMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
axis(4)
mtext("Raup-Crick distance", side=4, line=2.8)
par(new=T)
plot(xval, BinaryRaupCrickDateMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickDateMCA$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BinaryRaupCrickDateMCA")
abline(v=BinaryRaupCrickDateMCA$breaks, lty=2)
dev.off()
