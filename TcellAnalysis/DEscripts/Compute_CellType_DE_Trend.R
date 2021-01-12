library(edgeR)
library(scran)
library(scater)

# ---- Load Data ----

sce.full <- readRDS("../../Newcastle/Cambridge_Sanger_Newcastle-multiBatchnorm.RDS")
tcells <- read.table("../data/Tcell_annotations_ext.tsv",sep="\t",header=TRUE)
colnames(tcells)[2] <- "Annotation"

# I am removing IgH and IgL genes, these are most-likely background 
rmGenes <- grep("^IGH|^IGL|^IGK",rownames(sce.full),value=T)
sce.full <- sce.full[!rownames(sce.full) %in% rmGenes,]

# Doesn't do character subsetting?
sce.full <- sce.full[,colnames(sce.full) %in% tcells$CellID]

# Add meta data from Gary
meta <- read.table("/hps/research/sds/sds-marioni-cov/Newcastle/combined_dec_MetaData.txt",sep="\t",header=TRUE,stringsAsFactor=FALSE)
meta <- as.data.frame(apply(meta,2, function(COL) gsub("^b","",COL)))

# Remove all previous meta data because that was wrong
meta$CellID <- rownames(meta) <- meta$X
meta <- meta[colnames(sce.full),]
meta <- dplyr::left_join(meta,tcells)
colData(sce.full) <- DataFrame(meta)

# Rm patient with B-Cell Lymphoma
bc <- c("CV0198")
sce.full <- sce.full[, sce.full$sample_id !="BGCV01_CV0902" & sce.full$patient_id != bc & !(sce.full$Status_on_day_collection_summary %in% c("Non_covid","LPS")) & sce.full$Collection_Day =="D0"]

# Rm remaining doublets
mk <- read.table("../../Cydar/data/doublets/Tcell_doublets.tsv",header=FALSE,sep="\t",stringsAsFactors=FALSE)
mk2 <- read.table("../../Cydar/data/doublets/Tcell_doublets_ext.tsv",header=FALSE,sep="\t",stringsAsFactors=FALSE)
sce.full <- sce.full[,!sce.full$CellID %in% c(mk$V1,mk$V2)]

colData(sce.full) <- droplevels(colData(sce.full))


# ---- Subset to Group of Interest ----

# Starting with Asymptomatic vs. Healthy (this is only for historic reasons, if you want to exclude critical or something)
sce <- sce.full[,sce.full$Status_on_day_collection_summary %in% c("Healthy","Asymptomatic","Mild","Moderate","Severe","Critical")]

ctpp <- table(sce$Annotation,sce$patient_id)
meta <- data.frame(colData(sce)[,c("patient_id","Status_on_day_collection_summary")])
meta <- meta[!duplicated(meta$patient_id),]

# I am requiring to have at least 3 patients in each group
minPatPerGroup <- 3
isMoreThanMin <- apply(ctpp,1, function(RW) {
	  pats <- colnames(ctpp)[RW<20]
	  meta.after <- meta[!meta$patient_id %in% pats,]
	  nHel <- sum(meta.after$Status_on_day_collection_summary=="Healthy")
	  nAsy <- sum(meta.after$Status_on_day_collection_summary=="Asymptomatic")
	  nMil <- sum(meta.after$Status_on_day_collection_summary=="Mild")
	  nMod <- sum(meta.after$Status_on_day_collection_summary=="Moderate")
	  nSev <- sum(meta.after$Status_on_day_collection_summary=="Severe")
	  nCrit <- sum(meta.after$Status_on_day_collection_summary=="Critical")
	  isMinNum <- all(c(nHel,nAsy,nMil,nMod,nSev,nCrit) >= 3)
	  return(isMinNum)
})

keepCells <- names(isMoreThanMin[isMoreThanMin])

# Create final sce
sce <- sce[,sce$Annotation %in% keepCells]

# ---- Loop over Cell Types ----

cts <- unique(sce$Annotation)

cfs <- c("Severity.L","Severity.Q")
for (ct in cts) {
    # Set Cell Type
    sce.sub <- sce[,sce$Annotation==ct]

    # Remove patients with less than 20 cells
    nCells <- table(sce.sub$patient_id)
    rmPatients <- names(nCells[nCells<20])
    sce.sub <- sce.sub[,!sce.sub$patient_id %in% rmPatients]
    print(rmPatients)

    # Summarize Counts
    smrzd <- aggregateAcrossCells(sce.sub,
				  id=as.character(colData(sce.sub)[,c("patient_id")]))

    y <- DGEList(counts=counts(smrzd),
		 samples=colData(smrzd))
    keep <- filterByExpr(y, group=y$samples$Status_on_day_collection_summary, min.count=3, min.total.count=5)
    y <- y[keep,]

    # Norm
    y <- calcNormFactors(y)

    # Set reference level
    y$samples$Severity <- factor(y$samples$Status_on_day_collection_summary,levels=c("Healthy","Asymptomatic",
								      "Mild","Moderate",
								      "Severe","Critical"))
    y$samples$Severity <- ordered(y$samples$Severity)
    y$samples <- droplevels(y$samples)
    y$samples$Age <- as.numeric(y$samples$Age)
    mdl <- model.matrix(~ 0 + Site + Age + Sex + Severity, data=y$samples)

    # Estimate Dispersion
    y <- estimateDisp(y, mdl)

    # Fit Model
    fit <- glmQLFit(y, mdl, robust=TRUE)

    for (cf in cfs) {
	res <- glmQLFTest(fit, coef=cf)

	# Table
	out <- topTags(res,n=nrow(y))$table
	out$CellType <- ct
	out$Gene <- rownames(out)
	out$Comp <- cf
	write.csv(out,paste0("../data/DE_",ct,"_",cf,".csv"))
    }
}
