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
dir.create("OverlappedPairedEnd_wSTD_13_RAnalysisResults")

# Make species-level barplot
pdf("OverlappedPairedEnd_wSTD_13_RAnalysisResults/barplottop50species.pdf", width=14, height=10)
top50species <- read.delim("OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_top50species_nreads_fishes_concentration.tsv", header=T, check.names=F)
temp <- ggplot(top50species, aes(x=samplename, y=nreads, fill=fct_rev(species)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="species")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
temp <- temp + ylab("number of copies / mL")
plot(temp)
dev.off()

# Make family-level barplot
pdf("OverlappedPairedEnd_wSTD_13_RAnalysisResults/barplottop50family.pdf", width=13, height=10)
top50family <- read.delim("OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_top50family_nreads_fishes_concentration.tsv", header=T, check.names=F)
temp <- ggplot(top50family, aes(x=samplename, y=nreads, fill=fct_rev(family)))
temp <- temp + geom_bar(stat="identity", position="fill")
temp <- temp + scale_y_continuous(labels=percent)
temp <- temp + scale_fill_manual(values=c("#C0C0C0FF", pal_igv(alpha=0.8)(50)), name="family")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
temp <- temp + ylab("number of copies / mL")
plot(temp)
dev.off()

# Make species-level heatmap
pdf("OverlappedPairedEnd_wSTD_13_RAnalysisResults/heatmapspecies.pdf", width=45, height=10)
commspecies <- read.delim("OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_species_nreads_fishes_concentration.tsv", header=T, check.names=F)
commspecies$nreads[(commspecies$nreads == 0)] <- NA
temp <- ggplot(commspecies, aes(x=species, y=samplename, fill=nreads))
temp <- temp + geom_tile()
temp <- temp + scale_fill_gsea(na.value="white")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
temp <- temp + labs(fill="number of copies / mL")
plot(temp)
dev.off()

# Make family-level heatmap
pdf("OverlappedPairedEnd_wSTD_13_RAnalysisResults/heatmapfamily.pdf", width=30, height=10)
commfamily <- read.delim("OverlappedPairedEnd_wSTD_12_ClaidentResults/sample_family_nreads_fishes_concentration.tsv", header=T, check.names=F)
commfamily$nreads[(commfamily$nreads == 0)] <- NA
temp <- ggplot(commfamily, aes(x=family, y=samplename, fill=nreads))
temp <- temp + geom_tile()
temp <- temp + scale_fill_gsea(na.value="white")
temp <- temp + theme_test()
temp <- temp + theme(axis.text.x=element_text(angle=90, hjust=1))
temp <- temp + labs(fill="number of copies / mL")
plot(temp)
dev.off()
