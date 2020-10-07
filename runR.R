library(tidyverse)
library(ggsci)
library(vegan)
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
# Coverage-based rarefaction

# Spatial Mantel correlogram analysis

# Temporal Mantel correlogram analysis

# PERMANOVA

# NMDS

# Phylogenetic community ecological analyses

