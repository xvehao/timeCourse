rm(list=ls())
library(DESeq2)
library(tidyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(reshape2)
library(annotate)
source('eda/eda_function.R')
hours = c(0,3,6,12,24,48,72)

load('../data/pathway/KEGG_common_genes.RData')
H3K27ac <- read.delim("../data/count/H3K27ac_merged.txt")
atac <- read.delim("../data/count/ATACseq_merged.tab")
rna <- read.delim("../data/count/RNAseq_merged.txt")
rnaDE = apply_deseq2(rna)
H3K27acDE = apply_deseq2(H3K27ac)
atacDE = apply_deseq2(atac)
clu = read_clu()

chip_gr = makeGRangesFromDataFrame(H3K27acDE$lfc[,c('Chr','Start','End')])
chip_peakAnno_df = annotate_peak(chip_gr)
atac_gr = makeGRangesFromDataFrame(atacDE$lfc[,c('Chr','Start','End')])
atac_peakAnno_df = annotate_peak(atac_gr)
### subset chip-seq ###
H3K27acDE$lfc = H3K27acDE$lfc %>% left_join(chip_peakAnno_df[,c('seqnames','start','end','GENENAME')],
                                            c('Start'='start','End'='end','Chr'='seqnames'))
H3K27acDE$count_normalized = H3K27acDE$count_normalized %>% left_join(chip_peakAnno_df[,c('seqnames','start','end','GENENAME')],
                                                                      c('Start'='start','End'='end','Chr'='seqnames'))
selected_H3K27ac_region = H3K27acDE$lfc$Start %in% clu$chip$start
selected_gene = rnaDE$lfc$GENENAME %in% clu$rna$gene
chip_lfc_selected = H3K27acDE$lfc[selected_H3K27ac_region,]
chip_lfc_selected = chip_lfc_selected[complete.cases(chip_lfc_selected),]
chip_count_selected = H3K27acDE$count_normalized[selected_H3K27ac_region,]
chip_count_selected = chip_count_selected[complete.cases(chip_count_selected),]
### subset atac-seq ###
atacDE$lfc = atacDE$lfc %>% left_join(atac_peakAnno_df[,c('seqnames','start','end','GENENAME')],
                                      c('Start'='start','End'='end','Chr'='seqnames'))
atacDE$count_normalized = atacDE$count_normalized %>% left_join(atac_peakAnno_df[,c('seqnames','start','end','GENENAME')],
                                                                c('Start'='start','End'='end','Chr'='seqnames'))
selected_atac_region = atacDE$lfc$Start %in% clu$atac$start
atac_lfc_selected = atacDE$lfc[selected_atac_region,]
atac_lfc_selected = atac_lfc_selected[complete.cases(atac_lfc_selected),]
atac_count_selected = atacDE$count_normalized[selected_atac_region,]
atac_count_selected = atac_count_selected[complete.cases(atac_count_selected),]
### subset rna-seq ###
rna_lfc_selected = rnaDE$lfc[selected_gene,]
rna_lfc_selected = rna_lfc_selected[complete.cases(rna_lfc_selected),]
rna_count_selected = rnaDE$count_normalized[selected_gene,]
rna_count_selected = rna_count_selected[complete.cases(rna_count_selected),]

genes_rna_chip = intersect(rna_lfc_selected$GENENAME,chip_lfc_selected$GENENAME)
genes_chip_atac = intersect(atac_lfc_selected$GENENAME,chip_lfc_selected$GENENAME)
genes_rna_atac = intersect(rna_lfc_selected$GENENAME,atac_lfc_selected$GENENAME)
genes_intersection = intersect(intersect(rna_lfc_selected$GENENAME,
                                         atac_lfc_selected$GENENAME),
                               chip_lfc_selected$GENENAME)
all_genes = unique(c(rna_lfc_selected$GENENAME,
                     chip_lfc_selected$GENENAME,
                     atac_lfc_selected$GENENAME))



GENE.DESCRIPTION.HUMAN_0 <- read.delim("../data/GENE-DESCRIPTION-TXT_HUMAN_0.txt", header=FALSE, comment.char="#")
n_gene = nrow(GENE.DESCRIPTION.HUMAN_0)
genenames = GENE.DESCRIPTION.HUMAN_0[seq(1,n_gene,by=2),2]
genedesp = GENE.DESCRIPTION.HUMAN_0[seq(2,n_gene,by=2),1]
GENE.DESCRIPTION = data.frame(genedesp)
rownames(GENE.DESCRIPTION)=genenames

destination = './eda/rna_chip_atac_pathway.pdf'
pdf(file=destination, onefile=T)
par(mfrow = c(1,1))
for (j in 1:nrow(KEGG)){
  genes = getSYMBOL(str_split(KEGG$geneID[j],'/')[[1]], data='org.Hs.eg')
  if((sum(rnaDE$lfc$GENENAME%in%genes)==0)|(sum(H3K27acDE$lfc$GENENAME%in%genes)==0)|(sum(atacDE$lfc$GENENAME%in%genes)==0)){
    print('genes skipped')
    next
  }
  gFC = plot_lfc(genes,lfcdat=rnaDE$lfc,ylab="RNA lfc")
  chipFC = plot_lfc(genes,lfcdat=H3K27acDE$lfc,ylab="ChIP lfc")
  atacFC = plot_lfc(genes,lfcdat=atacDE$lfc,ylab="ATAC lfc")
  grid.arrange(gFC, chipFC, atacFC, nrow=2)
}
dev.off()

source("Bayesian-ODE-Models-Genetics/1_LLR2.R")
rownames(chip_lfc_selected) = paste0('chip_',make.names(chip_lfc_selected$GENENAME,unique=T))
geneData <- rbind(rna_lfc_selected,chip_lfc_selected[,1:7])
siteNames <- rownames(geneData)
# ========================================================================
nonBayesLLR2 <- LLR2(geneData, hours, bayes=FALSE, writeToCSV=TRUE)