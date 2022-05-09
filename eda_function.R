library(DESeq2)
library(genomation)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)

apply_deseq2 = function(dat, hours=c(0,3,6,12,24,48,72)){
  exp_level_col = str_detect(colnames(dat),pattern='hr_')
  n_rep = sum(exp_level_col)/length(hours)
  dat_rmNA = dat[complete.cases(dat),]
  coldat = data.frame(condition=factor(rep(c("untreated","treated"),each=n_rep)))
  coldat$condition <- relevel(coldat$condition, ref = "untreated")
  lfc = data.frame(matrix(0,nrow=nrow(dat_rmNA),ncol=length(hours)))
  pvalue = padj = data.frame(matrix(1,nrow=nrow(dat_rmNA),ncol=length(hours)))
  colnames(lfc) = colnames(pvalue) = colnames(padj)=c(paste0(hours,'hr'))
  count_normalized = dat_rmNA
  ref_col = colnames(dat_rmNA)[str_detect(colnames(dat_rmNA),pattern='0hr')]
  for (h in hours[-1]){
    trt_col = colnames(dat_rmNA)[str_detect(colnames(dat_rmNA),pattern=paste0(h,'hr'))]
    dds <- DESeqDataSetFromMatrix(countData = dat_rmNA[,c(ref_col,trt_col)],
                                  colData = coldat,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res = results(dds)
    lfc[paste0(h,'hr')]=res$log2FoldChange
    pvalue[paste0(h,'hr')]=res$pvalue
    padj[paste0(h,'hr')]=res$padj
    count_normalized[,trt_col] = counts(dds)[,trt_col]
  }
  count_normalized[,ref_col] = counts(dds)[,ref_col]
  if('geneName' %in% colnames(dat_rmNA)){
    lfc['GENENAME'] = pvalue['geneName'] = padj['geneName'] = dat_rmNA$geneName
  }else{
    loc_idx = colnames(dat_rmNA)[!exp_level_col]
    lfc[,loc_idx] = pvalue[,loc_idx] = padj[,loc_idx] = dat_rmNA[,loc_idx]
  }
  return(list(lfc=lfc,pvalue=pvalue,padj=padj,count_normalized=count_normalized))
}

annotate_peak = function(peak_gr){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  peakAnno <- annotatePeak(peak_gr, TxDb=txdb, verbose=FALSE)
  peakAnno_df = data.frame(peakAnno)
  entrez <- peakAnno_df$geneId
  
  # Return the gene symbol for the set of Entrez IDs
  annotations_edb <- AnnotationDbi::select(EnsDb.Hsapiens.v75,
                                           keys = entrez,
                                           columns = c("GENENAME"),
                                           keytype = "ENTREZID")
  annotations_edb$ENTREZID <- as.character(annotations_edb$ENTREZID)
  peakAnno_DF = peakAnno_df %>% 
    left_join(annotations_edb, by=c("geneId"="ENTREZID"))
  return(peakAnno_DF)
}

read_clu = function(){
  mmc3 <- read_excel("../references/atac-seq/mmc3.xlsx",
                     sheet = "2.Differential peaks-expression")
  chip = mmc3[-1,1:4]
  colnames(chip) = mmc3[1,1:4]
  atac = mmc3[-1,6:9]
  colnames(atac) = mmc3[1,6:9]
  rna = mmc3[-1,11:12]
  colnames(rna) = mmc3[1,11:12]
  return(list(chip=chip,rna=rna,atac=atac))
}

plot_lfc = function(genes,lfcdat=rnaDE$lfc,ylab="RNA lfc"){
  datFC = lfcdat[lfcdat$GENENAME %in% genes,]
  # fold change plot:
  datFC = gather (datFC, hours, FC, '0hr':'72hr', factor_key=TRUE)
  datFC$hours = gsub('hr', '', datFC$hours)
  if ('Start' %in% colnames(datFC)){
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           group=Start,
                           color=GENENAME
    ))+
      geom_line(size=.5, alpha=0.4) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      ggtitle(paste0(c(KEGG[j,c('ID','Description','GeneRatio','p.adjust','Count')]),collapse = ' '))+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }else{
    ggplot(data=datFC, aes(x=as.numeric(as.character(hours)),
                           y=as.numeric(FC), 
                           color=GENENAME
    ))+
      geom_line(size=.5) + geom_point(size=1) +
      theme_bw()+
      theme(strip.text.x = element_text(size=5),
            plot.title=element_text(size=5),
            axis.text=element_text(size=5),
            legend.text=element_text(size=5),
            legend.title=element_text(size=5),
            axis.title.x=element_text(size=5, margin=margin(15,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,5,0,0))) +
      scale_x_continuous("time point number", breaks=as.numeric((datFC$hours)), labels =
                           as.numeric(as.character(datFC$hours)),
                         sec.axis = sec_axis(~., breaks=as.numeric((datFC$hours)), labels =
                                               as.numeric(as.character(datFC$hours)),
                                             name ="hours"))+
      ylab(ylab)+
      ggtitle(paste0(c(KEGG[j,c('ID','Description','GeneRatio','p.adjust','Count')]),collapse = ' '))+
      theme(plot.title = element_textbox_simple(
        size = 5, lineheight = -1, padding = margin(0, 0, 5, 0)))
  }
}