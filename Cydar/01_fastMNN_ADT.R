library(batchelor)
library(scater)
library(BiocParallel)
library(BiocSingular)

adtf <- "../../Newcastle/combined_dec_ADT_SCE.RDS"

adt <- readRDS(adtf)

# There are some samples with no antibody signal, cool
rmSamples <- c("MH8919226","MH8919227","MH8919228","MH8919229","MH8919230","MH8919231","MH8919232","MH8919233")

# removing doublets identified by Mike and Kelvin
mk <- read.table("../data/doublets/Tcell_doublets.tsv",header=FALSE,sep="\t",stringsAsFactors=FALSE)
kl <- read.csv("../data/doublets/Karsten_doublets.csv",stringsAsFactors=FALSE)
dbls <- c(kl$doublets[kl$doublets!=""], mk$V1)

rmBcs <- colnames(adt)[adt$initial_clustering=="Doublet" | adt$sample_id %in% rmSamples]
rmBcs <- c(rmBcs,dbls)



adt <- adt[,!colnames(adt) %in% rmBcs]

set.seed(42)
adt1 <- adt[,adt$Site=="Cambridge"]
adt2 <- adt[,adt$Site=="Ncl"]
adt3 <- adt[,adt$Site=="Sanger"]


param <- MulticoreParam(workers=4)
print("Starting MNN")
set.seed(300)
mnncor <- batchelor::fastMNN(adt2,adt1,adt3,
		     BPPARAM=param, # commented this out for singularity
		     k=20,
		     d=50,
		     merge.order=c(1,2,3),
		     BSPARAM=IrlbaParam(deferred=TRUE),
		     assay.type="counts",
		     #                      BNPARAM=AnnoyParam(),
		     cos.norm=TRUE, # prob necesarry
		     correct.all=TRUE)
print("Done MNN")

m.cor <- assay(mnncor,"reconstructed")
m.cor <- m.cor[,colnames(adt)]
assay(adt,"reconstructed") <- m.cor
saveRDS(adt,"../data/combined_dec_ADT_MNNcorrected.RDS")
