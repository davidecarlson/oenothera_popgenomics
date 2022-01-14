#!/usr/bin/env Rscript

library(gdsfmt)
library(SNPRelate)
library(tidyverse)
library(pheatmap)

# take in vcf file and produce distance based dendrogram

vcf <- '../vcfs/variants.32samples.0.5missing.biallelicSNPs.vcf.gz'



# convert VCF to GDS format

snpgdsVCF2GDS(vcf, "../vcfs/data.gds")

# summarize the results

snpgdsSummary("../vcfs/data.gds")

genotypes<-snpgdsOpen("../vcfs/data.gds")
#pop_code <- scan("../popmap_nogrand.txt", what=character())

pop_code <- c("biennis","biennis","biennis","biennis","biennis","biennis","biennis",
            "biennis","biennis","biennis","biennis","biennis","biennis","biennis",
            "elata","elata","elata","elata","elata","elata","elata","elata","elata",
            "elata","elata","elata","elata","elata","elata","biennis","biennis","elata")
pops <- as.factor(pop_code)

print(pop_code)

# get distance matrix and dendogram

dissMatrix  <-  snpgdsDiss(genotypes , sample.id=NULL, snp.id=NULL, autosome.only=FALSE,
                remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, num.thread=10, verbose=TRUE)


# Hierarchical clustering
snpclust <- snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)

# make a dendrogram

snpTree <- snpgdsCutTree(snpclust, samp.group=as.factor(pop_code))
        
#print(snpTree)

pdf('../figures/snp_dendro.pdf')

plot(snpTree$dendrogram)
legend("topright", legend=levels(pops), col=1:nlevels(pops), pch=12, ncol=2)

dev.off()


# get IBS results and also use that to produce a dendrogram

ibs.hc <- snpgdsHCluster(snpgdsIBS(genotypes, num.thread=2, autosome.only = FALSE))

# individulas in the same population are clustered together





ibsTree <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

print(ibsTree)


pdf('../figures/ibs_dendro.pdf')

plot(ibsTree$dendrogram)
legend("topright", legend=levels(pops), col=1:nlevels(pops), pch=12, ncol=2)

dev.off()
