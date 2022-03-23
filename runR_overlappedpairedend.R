library(tidyverse)
library(ggsci)
library(foreach)
library(doParallel)
library(vegan)
library(mpmcorrelogram)
library(geosphere)
library(scales)
library(pvclust)

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
pdf("OverlappedPairedEnd_11_RAnalysisResults/heatmapspecies.pdf", width=45, height=10)
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
pdf("OverlappedPairedEnd_11_RAnalysisResults/heatmapfamily.pdf", width=30, height=10)
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
Community <- read.table("OverlappedPairedEnd_10_ClaidentResults/sample_otu_matrix_fishes.tsv", header=T, row.names=1)

# Draw OTU accumulation curve
SpecAccum <- specaccum(Community)
pdf("OverlappedPairedEnd_11_RAnalysisResults/specaccum.pdf", width=7, height=7)
plot(SpecAccum, xlab="number of samples", ylab="number of OTUs", main="OTU accumulation curve")
dev.off()

# Draw rarefaction curves
pdf("OverlappedPairedEnd_11_RAnalysisResults/rarecurve.pdf", width=7, height=7)
rarecurve(Community, step=10, xlab="number of seqs", ylab="number of OTUs", main="rarefaction curves")
dev.off()

# Coverage-based rarefaction
## make rareslopelist using all cpu cores
rareslopelist <- list()
cl <- makeCluster(detectCores(), type="FORK")
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
temp <- as.data.frame(row.names(Community), row.names=row.names(Community))
colnames(temp) <- "samplename"
write.table(cbind(temp, as.data.frame(cvrrare)), "OverlappedPairedEnd_11_RAnalysisResults/cvrrare.tsv", sep="\t", append=F, quote=F, row.names=F, col.names=F, na="NA")
RarefiedCommunity <- list()
for(i in 1:4) {
  RarefiedCommunity[[i]] <- rrarefy(Community, cvrrare)
  write.table(cbind(temp, RarefiedCommunity[[i]]), paste0("OverlappedPairedEnd_11_RAnalysisResults/RarefiedCommunity", i, ".tsv"), sep="\t", append=F, quote=F, row.names=F, col.names=T, na="NA")
}

# Make Beta-diversity (dissimilarity) matrix
BrayCurtis <- list()
Jaccard <- list()
BinaryJaccard <- list()
BinaryRaupCrick <- list()
for(i in 1:4) {
  BrayCurtis[[i]] <- vegdist(RarefiedCommunity[[i]], method="bray")
  Jaccard[[i]] <- vegdist(RarefiedCommunity[[i]], method="jaccard")
  BinaryJaccard[[i]] <- vegdist(RarefiedCommunity[[i]], method="jaccard", binary=T)
  BinaryRaupCrick[[i]] <- as.dist(raupcrick(RarefiedCommunity[[i]], null="r1", nsimul=999))
}

# Read metadata
Metadata <- read.table("Metadata.tsv", header=T, row.names=1)

# PERMANOVA
sink("OverlappedPairedEnd_11_RAnalysisResults/PERMANOVA.txt", split=T)
for(i in 1:4) {
  print(adonis(BrayCurtis[[i]] ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores()))
}
for(i in 1:4) {
  print(adonis(Jaccard[[i]] ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores()))
}
for(i in 1:4) {
  print(adonis(BinaryJaccard[[i]] ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores()))
}
for(i in 1:4) {
  print(adonis(BinaryRaupCrick[[i]] ~ as.factor(Metadata$Type) + as.numeric(Metadata$Temperature) + as.numeric(Metadata$Latitude) + as.factor(Metadata$Month) + 1, permutations=9999, parallel=detectCores()))
}
sink()

# Cluster analysis
BrayCurtisClusterSites <- list()
JaccardClusterSites <- list()
BinaryJaccardClusterSites <- list()
BinaryRaupCrickClusterSites <- list()
EuclideanClusterSpecies <- list()
BinaryEuclideanClusterSpecies <- list()
for(i in 1:4) {
  BrayCurtisClusterSites[[i]] <- pvclust(as.data.frame(t(RarefiedCommunity[[i]])), method.hclust="average", method.dist=function(x){vegan::vegdist(as.data.frame(t(x)),method="bray")}, nboot=100, parallel=T)
  JaccardClusterSites[[i]] <- pvclust(as.data.frame(t(RarefiedCommunity[[i]])), method.hclust="average", method.dist=function(x){vegan::vegdist(as.data.frame(t(x)),method="jaccard")}, nboot=100, parallel=T)
  BinaryJaccardClusterSites[[i]] <- pvclust(as.data.frame(t(RarefiedCommunity[[i]])), method.hclust="average", method.dist=function(x){vegan::vegdist(as.data.frame(t(x)),method="jaccard",binary=T)}, nboot=100, parallel=T)
  BinaryRaupCrickClusterSites[[i]] <- pvclust(as.data.frame(t(RarefiedCommunity[[i]])), method.hclust="average", method.dist=function(x){as.dist(vegan::raupcrick(as.data.frame(t(x)),null="r1",nsimul=999))}, nboot=100, parallel=T)
  EuclideanClusterSpecies[[i]] <- pvclust(as.data.frame(RarefiedCommunity[[i]]), method.hclust="average", method.dist=function(x){vegan::vegdist(as.data.frame(t(x)),method="euclidean")}, nboot=100, parallel=T)
  BinaryEuclideanClusterSpecies[[i]] <- pvclust(as.data.frame(RarefiedCommunity[[i]]), method.hclust="average", method.dist=function(x){vegan::vegdist(as.data.frame(t(x)),method="euclidean",binary=T)}, nboot=100, parallel=T)
}
## draw dendrograms
pdf("OverlappedPairedEnd_11_RAnalysisResults/ClusterAnalysis_sites.pdf", width=7, height=7)
for(i in 1:4) {
  plot(BrayCurtisClusterSites[[i]], xlab="site", ylab="Bray-Curtis distance", main="Cluster analysis among sites")
}
for(i in 1:4) {
  plot(JaccardClusterSites[[i]], xlab="site", ylab="Jaccard distance", main="Cluster analysis among sites")
}
for(i in 1:4) {
  plot(BinaryJaccardClusterSites[[i]], xlab="site", ylab="Jaccard distance of binary-transformed data", main="Cluster analysis among sites")
}
for(i in 1:4) {
  plot(BinaryRaupCrickClusterSites[[i]], xlab="site", ylab="Raup-Crick distance", main="Cluster analysis among sites")
}
dev.off()
pdf("OverlappedPairedEnd_11_RAnalysisResults/ClusterAnalysis_species.pdf", width=49, height=7)
for(i in 1:4) {
  plot(EuclideanClusterSpecies[[i]], xlab="site", ylab="Euclidean distance", main="Cluster analysis among species")
}
for(i in 1:4) {
  plot(BinaryEuclideanClusterSpecies[[i]], xlab="site", ylab="Euclidean distance of binary-transformed data", main="Cluster analysis among species")
}
dev.off()

# NMDS
## fit NMDS
BrayCurtisNMDS <- list()
JaccardNMDS <- list()
BinaryJaccardNMDS <- list()
BinaryRaupCrickNMDS <- list()
for(i in 1:4) {
  BrayCurtisNMDS[[i]] <- metaMDS(BrayCurtis[[i]], k=2, trymax=100)
  JaccardNMDS[[i]] <- metaMDS(Jaccard[[i]], k=2, trymax=100)
  BinaryJaccardNMDS[[i]] <- metaMDS(BinaryJaccard[[i]], k=2, trymax=100)
  BinaryRaupCrickNMDS[[i]] <- metaMDS(BinaryRaupCrick[[i]], k=2, trymax=100)
}
## renew NMDS by better-fitted results
for(i in 1:4) {
  for(j in 1:9) {
    set.seed(j)
    temp <- metaMDS(BrayCurtis[[i]], k=2, trymax=100)
    if(temp$stress < BrayCurtisNMDS[[i]]$stress) {
      BrayCurtisNMDS[[i]] <- temp
    }
    temp <- metaMDS(Jaccard[[i]], k=2, trymax=100)
    if(temp$stress < JaccardNMDS[[i]]$stress) {
      JaccardNMDS[[i]] <- temp
    }
    temp <- metaMDS(BinaryJaccard[[i]], k=2, trymax=100)
    if(temp$stress < BinaryJaccardNMDS[[i]]$stress) {
      BinaryJaccardNMDS[[i]] <- temp
    }
    temp <- metaMDS(BinaryRaupCrick[[i]], k=2, trymax=100)
    if(temp$stress < BinaryRaupCrickNMDS[[i]]$stress) {
      BinaryRaupCrickNMDS[[i]] <- temp
    }
  }
}
## fit environmental data
BrayCurtisNMDSenv <- list()
JaccardNMDSenv <- list()
BinaryJaccardNMDSenv <- list()
BinaryRaupCrickNMDSenv <- list()
for(i in 1:4) {
  BrayCurtisNMDSenv[[i]] <- envfit(BrayCurtisNMDS[[i]], Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
  JaccardNMDSenv[[i]] <- envfit(JaccardNMDS[[i]], Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
  BinaryJaccardNMDSenv[[i]] <- envfit(BinaryJaccardNMDS[[i]], Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
  BinaryRaupCrickNMDSenv[[i]] <- envfit(BinaryRaupCrickNMDS[[i]], Metadata[,c("Type", "Temperature", "Latitude", "Date")], permu=999)
}
## draw NMDS
pdf("OverlappedPairedEnd_11_RAnalysisResults/NMDS.pdf", width=7, height=7)
for(i in 1:4) {
  ordiplot(BrayCurtisNMDS[[i]], type="n")
  orditorp(BrayCurtisNMDS[[i]], display="sites", air=0.1, cex=1)
  plot(BrayCurtisNMDSenv[[i]], p.max=0.05)
}
for(i in 1:4) {
  ordiplot(JaccardNMDS[[i]], type="n")
  orditorp(JaccardNMDS[[i]], display="sites", air=0.1, cex=1)
  plot(JaccardNMDSenv[[i]], p.max=0.05)
}
for(i in 1:4) {
  ordiplot(BinaryJaccardNMDS[[i]], type="n")
  orditorp(BinaryJaccardNMDS[[i]], display="sites", air=0.1, cex=1)
  plot(BinaryJaccardNMDSenv[[i]], p.max=0.05)
}
for(i in 1:4) {
  ordiplot(BinaryRaupCrickNMDS[[i]], type="n")
  orditorp(BinaryRaupCrickNMDS[[i]], display="sites", air=0.1, cex=1)
  plot(BinaryRaupCrickNMDSenv[[i]], p.max=0.05)
}
dev.off()

# Spatial Mantel correlogram analysis
## make geographical distance matrix
Geodist <- as.dist(distm(cbind(Metadata$Longitude, Metadata$Latitude), fun=distGeo))
## run analysis
BrayCurtisGeoMCA <- list()
JaccardGeoMCA <- list()
BinaryJaccardGeoMCA <- list()
BinaryRaupCrickGeoMCA <- list()
for(i in 1:4) {
  BrayCurtisGeoMCA[[i]] <- mpmcorrelogram(BrayCurtis[[i]], Geodist, method="spearman", permutations=999, plot=F)
  JaccardGeoMCA[[i]] <- mpmcorrelogram(Jaccard[[i]], Geodist, method="spearman", permutations=999, plot=F)
  BinaryJaccardGeoMCA[[i]] <- mpmcorrelogram(BinaryJaccard[[i]], Geodist, method="spearman", permutations=999, plot=F)
  BinaryRaupCrickGeoMCA[[i]] <- mpmcorrelogram(BinaryRaupCrick[[i]], Geodist, method="spearman", permutations=999, plot=F)
}
## draw analysis results
pdf("OverlappedPairedEnd_11_RAnalysisResults/GeoMCA.pdf", width=7, height=7)
for(i in 1:4) {
  tipos <- BrayCurtisGeoMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BrayCurtisGeoMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BrayCurtisGeoMCA[[i]]$breaks[j] + BrayCurtisGeoMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Geodist, BrayCurtis[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BrayCurtisGeoMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Bray-Curtis distance", side=4, line=2.8)
  par(new=T)
  plot(xval, BrayCurtisGeoMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisGeoMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BrayCurtisGeoMCA")
  abline(v=BrayCurtisGeoMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- JaccardGeoMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(JaccardGeoMCA[[i]]$breaks) - 1)) {
    xval[j] <- (JaccardGeoMCA[[i]]$breaks[j] + JaccardGeoMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Geodist, Jaccard[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(JaccardGeoMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Jaccard distance", side=4, line=2.8)
  par(new=T)
  plot(xval, JaccardGeoMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardGeoMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="JaccardGeoMCA")
  abline(v=JaccardGeoMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- BinaryJaccardGeoMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BinaryJaccardGeoMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BinaryJaccardGeoMCA[[i]]$breaks[j] + BinaryJaccardGeoMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Geodist, BinaryJaccard[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BinaryJaccardGeoMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Jaccard distance of binary-transformed data", side=4, line=2.8)
  par(new=T)
  plot(xval, BinaryJaccardGeoMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardGeoMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BinaryJaccardGeoMCA")
  abline(v=BinaryJaccardGeoMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- BinaryRaupCrickGeoMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BinaryRaupCrickGeoMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BinaryRaupCrickGeoMCA[[i]]$breaks[j] + BinaryRaupCrickGeoMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Geodist, BinaryRaupCrick[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BinaryRaupCrickGeoMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Raup-Crick distance", side=4, line=2.8)
  par(new=T)
  plot(xval, BinaryRaupCrickGeoMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickGeoMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BinaryRaupCrickGeoMCA")
  abline(v=BinaryRaupCrickGeoMCA[[i]]$breaks, lty=2)
}
dev.off()

# Temporal Mantel correlogram analysis
## make date distance matrix
Datedist <- dist(as.Date(Metadata$Date))
## run analysis
BrayCurtisDateMCA <- list()
JaccardDateMCA <- list()
BinaryJaccardDateMCA <- list()
BinaryRaupCrickDateMCA <- list()
for(i in 1:4) {
  BrayCurtisDateMCA[[i]] <- mpmcorrelogram(BrayCurtis[[i]], Datedist, method="spearman", permutations=999, plot=F)
  JaccardDateMCA[[i]] <- mpmcorrelogram(Jaccard[[i]], Datedist, method="spearman", permutations=999, plot=F)
  BinaryJaccardDateMCA[[i]] <- mpmcorrelogram(BinaryJaccard[[i]], Datedist, method="spearman", permutations=999, plot=F)
  BinaryRaupCrickDateMCA[[i]] <- mpmcorrelogram(BinaryRaupCrick[[i]], Datedist, method="spearman", permutations=999, plot=F)
}
## draw analysis results
pdf("OverlappedPairedEnd_11_RAnalysisResults/DateMCA.pdf", width=7, height=7)
for(i in 1:4) {
  tipos <- BrayCurtisDateMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BrayCurtisDateMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BrayCurtisDateMCA[[i]]$breaks[j] + BrayCurtisDateMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Datedist, BrayCurtis[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BrayCurtisDateMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Bray-Curtis distance", side=4, line=2.8)
  par(new=T)
  plot(xval, BrayCurtisDateMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BrayCurtisDateMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BrayCurtisDateMCA")
  abline(v=BrayCurtisDateMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- JaccardDateMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(JaccardDateMCA[[i]]$breaks) - 1)) {
    xval[j] <- (JaccardDateMCA[[i]]$breaks[j] + JaccardDateMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Datedist, Jaccard[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(JaccardDateMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Jaccard distance", side=4, line=2.8)
  par(new=T)
  plot(xval, JaccardDateMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(JaccardDateMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="JaccardDateMCA")
  abline(v=JaccardDateMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- BinaryJaccardDateMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BinaryJaccardDateMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BinaryJaccardDateMCA[[i]]$breaks[j] + BinaryJaccardDateMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Datedist, BinaryJaccard[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BinaryJaccardDateMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Jaccard distance of binary-transformed data", side=4, line=2.8)
  par(new=T)
  plot(xval, BinaryJaccardDateMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryJaccardDateMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BinaryJaccardDateMCA")
  abline(v=BinaryJaccardDateMCA[[i]]$breaks, lty=2)
}
for(i in 1:4) {
  tipos <- BinaryRaupCrickDateMCA[[i]]$pval.Bonferroni < 0.05
  tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
  xval <- c()
  for(j in 1:(length(BinaryRaupCrickDateMCA[[i]]$breaks) - 1)) {
    xval[j] <- (BinaryRaupCrickDateMCA[[i]]$breaks[j] + BinaryRaupCrickDateMCA[[i]]$breaks[j + 1]) / 2
  }
  par(mar=c(5, 4, 4, 4))
  plot(Datedist, BinaryRaupCrick[[i]], pch=16, cex=1, type="p", xlim=c(0, tail(BinaryRaupCrickDateMCA[[i]]$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
  axis(4)
  mtext("Raup-Crick distance", side=4, line=2.8)
  par(new=T)
  plot(xval, BinaryRaupCrickDateMCA[[i]]$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(BinaryRaupCrickDateMCA[[i]]$breaks, n=1)), ylim=c(-1,1), xlab="date interval", ylab="Mantel correlation", main="BinaryRaupCrickDateMCA")
  abline(v=BinaryRaupCrickDateMCA[[i]]$breaks, lty=2)
}
dev.off()
