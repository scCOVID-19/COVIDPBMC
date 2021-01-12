#! /usr/bin/env Rscript

library(tidytext)
library(reshape2)
library(dplyr)
library(igraph)
library(stringdist)

message("Reading in clinical meta data")
covid.meta <- read.csv("/hps/research/sds/sds-marioni-cov/Newcastle/total_integrated.csv",
                       header=TRUE, stringsAsFactors=FALSE)

message("Reading in single-cell meta data")
all.meta <- read.table("/hps/research/sds/sds-marioni-cov/Newcastle/COVID19_scMeta-data.tsv",
                       sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(all.meta$CellID)

all.meta$Centre <- "Newcastle"
all.meta$Centre[grepl(all.meta$CellID, pattern="BGCV")] <- "Cambridge"
all.meta$Centre[grepl(all.meta$patient_id, pattern="AP")] <- "Sanger"

message("Reading in Sanger TCR data")
tcr.sanger <- read.table("/hps/research/sds/sds-marioni-cov/TCR_clonotype/Sanger_VDJ_merged.tsv",
                         sep="\t", header=TRUE, stringsAsFactors=FALSE)
tcr.sanger$Centre <- "Sanger"

message(paste0(nrow(tcr.sanger), " TCRs from Sanger"))

message("Reading in Newcastle TCR data")
tcr.ncl <- read.table("/hps/research/sds/sds-marioni-cov/TCR_clonotype/Newcastle_VDJ_merged.tsv",
                      sep="\t", header=TRUE, stringsAsFactors=FALSE)
tcr.ncl$Centre <- "Newcastle"

message(paste0(nrow(tcr.ncl), " TCRs from Newcastle"))

message("Reading in Cambridge TCR data")
tcr.cam <- read.table("//hps/research/sds/sds-marioni-cov/TCR_clonotype/Cambridge_VDJ_merged.tsv",
                      sep="\t", header=TRUE, stringsAsFactors=FALSE)
tcr.cam$Centre <- "Cambridge"

message(nrow(tcr.cam), " TCRs from Cambridge")

common.cols <- intersect(colnames(tcr.cam), intersect(colnames(tcr.sanger), colnames(tcr.ncl)))
tcr.all <- do.call(rbind.data.frame, list(tcr.sanger[, common.cols], 
                                          tcr.ncl[, common.cols], 
                                          tcr.cam[, common.cols]))

tcr.merge <- merge(tcr.all, all.meta, by=c('CellID', 'Centre'))
print(table(tcr.merge$Centre))

message("Saving full TCR data")
write.table(tcr.merge,
	    file="/hps/research/sds/sds-marioni-cov/TCR_clonotype/TCR_merged.tsv",
	    sep="\t", quote=FALSE, row.names=FALSE)

d.types <- unique(tcr.merge$patient_id)
chain.types <- unique(tcr.merge$chain)
dist.thresh <- 0

tcr.dist.list <- list()

message("Looping over all clonotypes based on Vgene + Jgene + CDR3aa length")
print(head(tcr.merge))
print(dim(tcr.merge))
tcr.merge$Text.Clonotype <- paste(tcr.merge$v_gene, tcr.merge$j_gene, nchar(tcr.merge$cdr3_nt), sep="_")

print(table(tcr.merge$Text.Clonotype, tcr.merge$Centre))

tcr.tab <- table(tcr.merge$Text.Clonotype)
mono.tcrs <- names(tcr.tab)[tcr.tab == 1]

mono.tcr.df <- tcr.merge[tcr.merge$Text.Clonotype %in% mono.tcrs, c("CellID", "Text.Clonotype")]
mono.tcr.df$TCR.Group <- paste(mono.tcr.df$Text.Clonotype, "G0", sep="_")

message(paste0("Found ", length(mono.tcrs), " singleton TCR clonotypes"))
t.types <- setdiff(unique(tcr.merge$Text.Clonotype), mono.tcrs)

message(paste0("Sub-grouping ", length(t.types), " non-singleton TCR clones"))
for(i in seq_along(t.types)){
    i.type <- t.types[i]
    sub.df <- tcr.merge[tcr.merge$Text.Clonotype %in% i.type, ]
    message(paste0("Clonotype ", i, ":", nrow(sub.df), " instances of ", i.type))
    if(nrow(sub.df) > 0){
        # within each group of identical text-based sequences, check for amino acid similarity (Hamming = 3).
        txt.clon.tab <- table(sub.df$Text.Clonotype)
        n.clons <- length(txt.clon.tab)
        multi.clons <- names(txt.clon.tab)[txt.clon.tab > 1]
        if(length(setdiff(names(txt.clon.tab), multi.clons))){
            uni.clons <- paste(setdiff(names(txt.clon.tab), multi.clons), "G0", sep="_")
            sub.groups <- c(uni.clons)
        } else{
            sub.groups <- c()
        }
        for(q in seq_along(multi.clons)){
            q.clon <- multi.clons[q]
            q.sub.dist <- (as.matrix(stringdistmatrix(a=sub.df$cdr3_nt[sub.df$Text.Clonotype %in% q.clon], method="hamming")) <= dist.thresh) + 0
            q.sub.graph <- graph_from_adjacency_matrix(q.sub.dist, mode="undirected")
            q.sub.groups <-  components(q.sub.graph)$membership
            # concatenate the group onto the clonotype
            q.sub.groups <- paste(sub.df$Text.Clonotype[sub.df$Text.Clonotype %in% q.clon], paste0("G", q.sub.groups), sep="_")
            sub.groups <- c(sub.groups, q.sub.groups)
        }
            
        tcr.dist.list[[i.type]] <- data.frame("CellID"=sub.df$CellID,
                                              "TCR.Group"=sub.groups,
                                              "Text.Clonotype"=i.type)
        sink(file="/dev/null")
        rm(list=c("sub.dist", "sub.groups", "sub.graph"))
        gc()
        sink(file=NULL)
    }
}

tcr.dist.df <- do.call(rbind.data.frame, tcr.dist.list)
tcr.dist.df <- do.call(rbind.data.frame, list(tcr.dist.df, mono.tcr.df[, colnames(tcr.dist.df)]))

message("Saving TCR clonotypes")
write.table(tcr.dist.df,
            file="/hps/research/sds/sds-marioni-cov/TCR_clonotype/TCR_clonotypes.tsv",
	    sep="\t", quote=FALSE, row.names=FALSE)