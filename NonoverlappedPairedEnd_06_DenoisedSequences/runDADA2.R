ranseed <- 1617568643
numthreads <- 32
outputfolder <- "NonoverlappedPairedEnd_06_DenoisedSequences"
pooling <- T
library(dada2)
set.seed(ranseed)
setDadaOpt(OMEGA_C=0)
fn1 <- sort(list.files(outputfolder, pattern="\\.fastq$", full.names=T))
extract.sample.names <- function (x) { sub("\\.fastq$", "", sub("^.*\\/", "", x)) }
names(fn1) <- sapply(fn1, extract.sample.names)
derep1 <- derepFastq(fn1, verbose=T, qualityType="FastqQuality")
err1 <- learnErrors(derep1, verbose=T, multithread=numthreads, qualityType="FastqQuality")
pdf(file=paste0(outputfolder, "/plotErrors.pdf"))
plotErrors(err1, obs=T, err_out=T, err_in=T, nominalQ=T)
dev.off()
dada1 <- dada(derep1, err=err1, verbose=T, multithread=numthreads, pool=pooling)
for (i in 1:length(fn1)) {
    write.table(derep1[[i]]$map, paste0(outputfolder, "/", names(fn1)[i], "_derepmap.txt"), sep="\t", col.names=F, row.names=F, quote=F)
    write.table(derep1[[i]]$uniques, paste0(outputfolder, "/", names(fn1)[i], "_uniques.txt"), sep="\t", col.names=F, row.names=T, quote=F)
    write.table(dada1[[i]]$map, paste0(outputfolder, "/", names(fn1)[i], "_dadamap.txt"), sep="\t", col.names=F, row.names=F, quote=F)
    write.table(dada1[[i]]$denoised, paste0(outputfolder, "/", names(fn1)[i], "_denoised.txt"), sep="\t", col.names=F, row.names=T, quote=F)
}
