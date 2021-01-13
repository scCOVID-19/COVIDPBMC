# ---- Load Data ----
library(scran)
library(scater)
library(randomForest)
library(ggplot2)
library(cowplot)
adtf <- "../data/combined_dec_ADT_MNNcorrected.RDS"
adt <- readRDS(adtf)

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-d", "--ds"), type="numeric",
       	             help="Fraction of cells to use for initating spheres")

opt <- parse_args(parser)

# ---- Subset to D0 cells ----

# Rm patient with B-Cell Lymphoma
bc <- c("CV0198")
adt <- adt[,adt$patient_id != bc & ! (adt$Status_on_day_collection_summary %in% c("Non_covid","LPS")) & adt$Collection_Day =="D0"]


# removing anything with less than 500 cells
drpSmps <- names(which(table(adt$patient_id)<500))
adt <- adt[,!adt$patient_id %in% drpSmps]

colData(adt) <- droplevels(colData(adt))

m.adt <- assay(adt,"reconstructed")

# ---- CYDAR ----
library(cydar)

# Create list
smps <- as.character(unique(adt$patient_id))


mrks <- rownames(adt)
mrks <- mrks[!grepl("Iso",mrks)]


mtx.list <- lapply(smps, function(SMP) {
	   cells <- colnames(adt)[adt$patient_id==SMP]
	   message(SMP)
	   m  <- t(m.adt[,cells])
	   return(m)
})
names(mtx.list) <- smps

cd <- prepareCellData(mtx.list,markers=mrks)

distances <- neighborDistances(cd,as.tol=TRUE)
req <- 20
med <- round(median(distances[,req-1]),2)

library(BiocParallel)
param <- MulticoreParam(workers=3)

set.seed(42)
n <- opt$ds
cd <- countCells(cd, tol=med, downsample=n, BPPARAM=param)

saveRDS(cd,paste0("../data/HyperSpheres_",n,".RDS"))
