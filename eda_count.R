rm(list=ls())
library(DESeq2)
library(tidyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(reshape2)
source('eda/eda_fun.R')
hours = c(0,3,6,12,24,48,72)

H3K27ac <- read.delim("../data/count/H3K27ac_merged.txt")
atac <- read.delim("../data/count/ATACseq_merged.tab")
rna <- read.delim("../data/count/RNAseq_merged.txt")
rnaDE = apply_deseq2(rna, n_rep=3)
H3K27acDE = apply_deseq2(H3K27ac, n_rep=2)
atacDE = apply_deseq2(atac, n_rep=2)
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
selected_gene = rownames(rnaDE$lfc) %in% clu$rna$gene
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

genes_rna_chip = intersect(rownames(rna_lfc_selected),chip_lfc_selected$GENENAME)
genes_chip_atac = intersect(atac_lfc_selected$GENENAME,chip_lfc_selected$GENENAME)
genes_rna_atac = intersect(rownames(rna_lfc_selected),atac_lfc_selected$GENENAME)
genes_intersection = intersect(intersect(rownames(rna_lfc_selected),
                                         atac_lfc_selected$GENENAME),
                               chip_lfc_selected$GENENAME)
all_genes = unique(c(rownames(rna_lfc_selected),
                     chip_lfc_selected$GENENAME,
                     atac_lfc_selected$GENENAME))



GENE.DESCRIPTION.HUMAN_0 <- read.delim("../data/GENE-DESCRIPTION-TXT_HUMAN_0.txt", header=FALSE, comment.char="#")
n_gene = nrow(GENE.DESCRIPTION.HUMAN_0)
genenames = GENE.DESCRIPTION.HUMAN_0[seq(1,n_gene,by=2),2]
genedesp = GENE.DESCRIPTION.HUMAN_0[seq(2,n_gene,by=2),1]
GENE.DESCRIPTION = data.frame(genedesp)
rownames(GENE.DESCRIPTION)=genenames
destination = './eda/rna_chip_atac.pdf'
pdf(file=destination, onefile=T)
par(mfrow = c(1,1))
for (j in 1:length(genes_intersection)){
  gene = genes_intersection[j]
  datFC = rna_lfc_selected[gene,]

  # fold change plot:
  datFC = gather (datFC, timepoint, FC, lfc0:lfc6, factor_key=TRUE)
  datFC$timepoint = 1:length(hours)
  datFC$hours = hours
  datFC$stage = c(rep("early", 5), rep("late", 2))
  genedescription = strsplit(GENE.DESCRIPTION[gene,],' ')[[1]]
  firsthalf_genedescription = paste0(genedescription[1:round(length(genedescription)/2)],collapse = ' ')
  secondhalf_genedescription = paste0(genedescription[(round(length(genedescription)/2)+1):length(genedescription)],collapse = ' ')
  gFC = ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                               y=as.numeric(FC)
  ))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                         as.numeric(as.character(datFC$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                             as.numeric(as.character(datFC$hours)),
                                           name ="hours"))+
    ylab("RNA log fold change")+
    ggtitle(firsthalf_genedescription)+
    theme(legend.position="none",
          plot.title = element_textbox_simple(
            size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))

  chipFC = chip_lfc_selected[chip_lfc_selected$GENENAME==gene,]
  n_peak = nrow(chipFC)
  chipFC = gather (chipFC, timepoint, FC, lfc0:lfc6, factor_key=TRUE)
  # fold change plot:
  chipFC$timepoint = rep(1:length(hours),each=n_peak)
  chipFC$hours = rep(hours,each=n_peak)
  chipFC$peaknum = rep(1:n_peak, length(hours))

  gFC_chip = ggplot(data=chipFC, aes(x=as.numeric(as.character(hours)),
                                         y=as.numeric(FC), colour = factor(peaknum),
  ))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                         as.numeric(as.character(datFC$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                             as.numeric(as.character(datFC$hours)),
                                           name ="hours"))+
    ylab("ChIP log fold change")+
    ggtitle(secondhalf_genedescription)+
    theme(legend.position="none",
          plot.title = element_textbox_simple(
            size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))

  # fold change plot:
  atacFC = atac_lfc_selected[atac_lfc_selected$GENENAME==gene,]
  n_peak = nrow(atacFC)
  atacFC = gather (atacFC, timepoint, FC, lfc0:lfc6, factor_key=TRUE)
  # fold change plot:
  atacFC$timepoint = rep(1:length(hours),each=n_peak)
  atacFC$hours = rep(hours,each=n_peak)
  atacFC$peaknum = rep(1:n_peak, length(hours))
  
  gFC_atac = ggplot(data=atacFC, aes(x=as.numeric(as.character(hours)),
                                     y=as.numeric(FC), colour = factor(peaknum),
  ))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((atacFC$hours)), labels =
                         as.numeric(as.character(atacFC$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric((atacFC$hours)), labels =
                                             as.numeric(as.character(atacFC$hours)),
                                           name ="hours"))+
    ylab("ATAC log fold change")+
    theme(legend.position="none",
          plot.title = element_textbox_simple(
            size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  
  # counts plot:
  datCounts = rna_count_selected[rna_count_selected$geneName==gene,]
  #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
  datCounts = gather (datCounts, timepoint, Counts, RNAseq_TP0hr_rep1:RNAseq_TP72hr_rep3, factor_key=TRUE)
  datCounts$timepoint = rep(c(1:length(hours)), each=3)
  datCounts$hours = rep(hours, each=3)
  datCounts$replicate = rep(c("A","B","C"),nrow(datCounts)/3)

  gCounts = ggplot(data=datCounts, aes(x=as.numeric(as.character(hours)),
                                       y=as.numeric(Counts), shape=replicate))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                         as.numeric(as.character(datFC$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                             as.numeric(as.character(datFC$hours)),
                                           name ="hours"))+
    ylab("RNA normalized counts")+
    ggtitle(paste0("normalized counts for both replicates ", gene))

  chipCounts = chip_count_selected[chip_count_selected$GENENAME==gene,]
  n_peak = nrow(chipCounts)
  #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
  chipCounts = gather (chipCounts, timepoint, Counts, K27ac_0hr_rep1:K27ac_72hr_rep2, factor_key=TRUE)
  chipCounts$timepoint = rep(c(1:length(hours)), each=2)
  chipCounts$hours = rep(hours, each=2*n_peak)
  chipCounts$replicate = factor(rep(c("A","B"),each=n_peak))
  chipCounts$peaknum = factor(rep(1:n_peak, length(hours)))
  gchipCounts = ggplot(data=chipCounts, aes(x=as.numeric(as.character(hours)),
                                            y=as.numeric(Counts),
                                            color=peaknum, shape=replicate))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((chipCounts$hours)), labels =
                         as.numeric(as.character(chipCounts$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric(hours), labels =
                                             as.numeric(as.character(hours)),
                                           name ="hours"))+
    ylab("ChIP normalized counts")+
    ggtitle(paste0("normalized counts for both replicates",gene))
  
  
  atacCounts = atac_count_selected[atac_count_selected$GENENAME==gene,]
  n_peak = nrow(atacCounts)
  #colnames(datCounts) = c("gene_ID", "gene_name", paste("X", colnames(datCounts)[-c(1,2)], sep=""))
  atacCounts = gather (atacCounts, timepoint, Counts, ATACseq_0hr_rep1:ATACseq_72hr_rep2, factor_key=TRUE)
  atacCounts$timepoint = rep(c(1:length(hours)), each=2)
  atacCounts$hours = rep(hours, each=2*n_peak)
  atacCounts$replicate = factor(rep(c("A","B"),each=n_peak))
  atacCounts$peaknum = factor(rep(1:n_peak, length(hours)))
  gatacCounts = ggplot(data=atacCounts, aes(x=as.numeric(as.character(hours)),
                                            y=as.numeric(Counts),
                                            color=peaknum, shape=replicate))+
    geom_line(size=1) + geom_point(size=3) +
    theme_bw()+
    theme(strip.text.x = element_text(size=5),
          plot.title=element_text(size=5),
          axis.text=element_text(size=5),
          legend.text=element_text(size=5),
          legend.title=element_text(size=5),
          axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
          axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
    scale_x_continuous("time point number", breaks=as.numeric((atacCounts$hours)), labels =
                         as.numeric(as.character(atacCounts$timepoint)),
                       sec.axis = sec_axis(~., breaks=as.numeric(hours), labels =
                                             as.numeric(as.character(hours)),
                                           name ="hours"))+
    ylab("atac normalized counts")+
    ggtitle(paste0("normalized counts for both replicates",gene))

  grid.arrange(gCounts, gchipCounts, gatacCounts, gFC, gFC_chip, gFC_atac, nrow=2)

}
dev.off()

source("Bayesian-ODE-Models-Genetics/1_LLR2.R")
rownames(chip_lfc_selected) = paste0('chip_',make.names(chip_lfc_selected$GENENAME,unique=T))
geneData <- rbind(rna_lfc_selected,chip_lfc_selected[,1:7])
siteNames <- rownames(geneData)
# ========================================================================
nonBayesLLR2 <- LLR2(geneData, hours, bayes=FALSE, writeToCSV=TRUE)