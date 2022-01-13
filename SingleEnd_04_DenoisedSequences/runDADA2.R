ranseed <- 1640907737
numthreads <- 128
outputfolder <- "SingleEnd_04_DenoisedSequences"
pooling <- T
library(dada2)
library(foreach)
library(doParallel)
set.seed(ranseed)
setDadaOpt(OMEGA_C=0)
fn1 <- sort(list.files(outputfolder, pattern="\\.fastq$", full.names=T))
extract.sample.names <- function (x) { sub("\\.fastq$", "", sub("^.*\\/", "", x)) }
names(fn1) <- sapply(fn1, extract.sample.names)
derep1 <- list()
cl <- makeCluster(numthreads, type="FORK")
registerDoParallel(cl)
derep1 <- foreach(i = 1:length(fn1), .packages="dada2") %dopar% {
    derepFastq(fn1[[i]], verbose=T, qualityType="FastqQuality")
}
stopCluster(cl)
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
