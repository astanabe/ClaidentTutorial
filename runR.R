library(tidyverse)
library(ggsci)
library(foreach)
library(doParallel)
library(vegan)
library(mpmcorrelogram)
library(geosphere)

# Make output directory
dir.create("11_RAnalysisResults")

# Make species-level Barplot
top10species <- read.table("08_FinalResults/sample_10species_matrix_bonyfishes.txt", header=T)
top10species$species <- factor(top10species$species, levels = c())
by(top10species, top10species$species, sum)
temp <- ggplot(top10species, aes(x=samplename, y=nreads, fill=species))
temp <- temp + geom_bar(stat = "identity", position = "fill")
temp <- temp + scale_y_continuous(labels = percent)
temp <- temp + scale_fill_manual(values=c(pal_igv()(10), "#C0C0C0FF"))
temp <- temp + theme_test()
temp <- temp + theme(axis.text = element_text(angle = 90))
plot(temp)

# Make family-level Barplot
top10family <- read.table("08_FinalResults/sample_10family_matrix_bonyfishes.txt", header=T)
temp <- ggplot(top10family, aes(x=samplename, y=nreads, fill=family))
temp <- temp + geom_bar(stat = "identity", position = "fill")
temp <- temp + scale_y_continuous(labels = percent)
temp <- temp + scale_fill_manual(values=c(pal_igv()(10), "#C0C0C0FF"))
temp <- temp + theme_test()
temp <- temp + theme(axis.text = element_text(angle = 90))
plot(temp)

# Read community data matrix
Community <- read.table("10_FinalResults/sample_species_matrix_fishes.tsv", header=T, row.names=1)

# Read metadata
Metadata <- read.table("Metadata.tsv", header=T, row.names=1)

# Draw species accumulation curve
SpecAccum <- specaccum(Community)
pdf("11_RAnalysisResults/specaccum.pdf")
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
#cvr <- max(getmincov)
cvr <- 0.01
## define function
cvrfun <- function(x) {min(which(x <= cvr))}
## get number of seqs of target coverage
cvrrare <- unlist(lapply(rareslopelist, cvrfun))
# make rarefied community data
RarefiedCommunity <- rrarefy(Community, cvrrare)
write.table(RarefiedCommunity, "11_RAnalysisResults/RarefiedCommunity.txt", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

# Make binary community data
BinaryRarefiedCommunity <- data.frame()
BinaryRarefiedCommunity <- replace(RarefiedCommunity, RarefiedCommunity > 0, 1)
write.table(BinaryRarefiedCommunity, "11_RAnalysisResults/BinaryRarefiedCommunity.txt", sep="\t", append=F, quote=F, row.names=T, col.names=T, na="NA")

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
pdf("11_RAnalysisResults/GeoMCA.pdf")
tipos <- BrayCurtisGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BrayCurtisGeoMCA$breaks[-1]), BrayCurtisGeoMCA$rM, pch=tipos, cex=1, type="b", xlab="geographic distance", ylab="Mantel correlation", main="BrayCurtisGeoMCA")
tipos <- JaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((JaccardGeoMCA$breaks[-1]), JaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlab="geographic distance", ylab="Mantel correlation", main="JaccardGeoMCA")
tipos <- BinaryJaccardGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BinaryJaccardGeoMCA$breaks[-1]), BinaryJaccardGeoMCA$rM, pch=tipos, cex=1, type="b", xlab="geographic distance", ylab="Mantel correlation", main="BinaryJaccardGeoMCA")
tipos <- BinaryRaupCrickGeoMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BinaryRaupCrickGeoMCA$breaks[-1]), BinaryRaupCrickGeoMCA$rM, pch=tipos, cex=1, type="b", xlab="geographic distance", ylab="Mantel correlation", main="BinaryRaupCrickGeoMCA")
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
pdf("11_RAnalysisResults/DateMCA.pdf")
tipos <- BrayCurtisDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BrayCurtisDateMCA$breaks[-1]), BrayCurtisDateMCA$rM, pch=tipos, cex=1, type="b", xlab="date interval", ylab="Mantel correlation", main="BrayCurtisDateMCA")
tipos <- JaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((JaccardDateMCA$breaks[-1]), JaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlab="date interval", ylab="Mantel correlation", main="JaccardDateMCA")
tipos <- BinaryJaccardDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BinaryJaccardDateMCA$breaks[-1]), BinaryJaccardDateMCA$rM, pch=tipos, cex=1, type="b", xlab="date interval", ylab="Mantel correlation", main="BinaryJaccardDateMCA")
tipos <- BinaryRaupCrickDateMCA$pval.Bonferroni < 0.05
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
plot((BinaryRaupCrickDateMCA$breaks[-1]), BinaryRaupCrickDateMCA$rM, pch=tipos, cex=1, type="b", xlab="date interval", ylab="Mantel correlation", main="BinaryRaupCrickDateMCA")
dev.off()
