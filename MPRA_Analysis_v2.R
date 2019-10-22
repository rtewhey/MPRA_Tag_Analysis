rm(list = ls())
library(DESeq2)
library(ggplot2)
options(stringsAsFactors = FALSE)

#############
###### CHANGE FOR EACH RUN
#############
wd<-"/home/user/wd/"

file_prefix="Project_Name"
file_date="20191002"
setwd(paste0(wd,"Analysis_dir/"))
dir.create("plots")
dir.create("results")

inputFullAttributesFileName=paste0(wd,"example.attributes")

countsFileName=paste0(wd,"Analysis_dir/example_counts.out")
inputConditionsFileName=paste0(wd,"Analysis_dir/example.cond")

#######

tag_counts<-read.table(countsFileName,header=TRUE)
fullattributeData=read.table(inputFullAttributesFileName, header=TRUE, row.names= NULL, check.names=FALSE,skipNul=TRUE)

if(!("haplotype" %in% colnames(fullattributeData)))
{
  fullattributeData$haplotype<-"NA"
}

counts<-read.table(file=countsFileName,header=TRUE)

#counts_oligo<-aggregate(. ~Oligo, data=tag_counts[,c(-1,-3,-4,-5,-6)], FUN=sum)
counts_oligo<-aggregate(. ~Oligo, data=tag_counts[,c(-1)], FUN=sum)

countData<-counts_oligo[,-1]
rownames(countData)<-counts_oligo[,1]
countDataRowNames=rownames(countData)

colData=read.table(inputConditionsFileName, header=FALSE, row.names=1)


colnames(colData)[1]="condition"
colData[,1]=factor(colData[,1])

colData$condition=relevel(colData$condition, "DNA")
rowNames=rownames(countData)
colnames(countData)<-row.names(colData)

#Reorder to fix prior cell line issue. Doing this incase any later code refers directly to column position
#Just doing this so the DF can be reordered....
colData$DNA<-0
colData[colData$condition=="DNA",]$DNA<-1 

################# initial processing of files #############



colData$condition=relevel(colData$condition, "DNA")
rowNames=rownames(countData)
colnames(countData)<-row.names(colData)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds$condition<-relevel(dds$condition, "DNA")
dds_results <- DESeq(dds,fitType='local')

orig_dds_results<-dds_results

#############
###### Normalization
#############

dds_results<-orig_dds_results

###Code for neg (all or ctrl) normalization
#correct on neg control
#nonexpA<-fullattributeData[fullattributeData$project=="negCtrl",]$ID

#correct from all
nonexpA<-fullattributeData$ID

#for (celltype in c("K562") {
for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  message(celltype)
  output_tmpA = results(dds_results, contrast=c("condition",celltype,"DNA"))
  
  ###Code for neg (all or ctrl) normalization
  #output_tmpA_neg<-output_tmpA[nonexpA,]
  #nonexpA<-row.names(output_tmpA_neg[!is.na(output_tmpA_neg$pvalue) & output_tmpA_neg$pvalue>0.001,])
  
  ###Code for summit shift normalization
  summit<-which.max(density(output_tmpA$log2FoldChange)$y)
  log_offset<-2^(density(output_tmpA$log2FoldChange)$x[summit])
  sizeFactors(dds_results)[which(colData$condition==celltype)]<-sizeFactors(dds_results)[which(colData$condition==celltype)]*(log_offset)
}

###Code for neg (all or ctrl) normalization
#dds_results_tmp<-estimateSizeFactors(dds[nonexpA])
#sizeFactors(dds_results)<-sizeFactors(dds_results_tmp)

###Code for summit shift normalization
dds_results<-estimateDispersions(dds_results,fitType='local')
dds_results<-nbinomWaldTest(dds_results)

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  
  output_tmpA = results(orig_dds_results, contrast=c("condition",celltype,"DNA"))
  outputA = results(dds_results, contrast=c("condition",celltype,"DNA"))
  
  pdf(paste0("plots/Normalized_FC_Density_",celltype,".pdf"),width=10,height=10)
  plot(density(output_tmpA[nonexpA,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="grey",main=paste0("Normalization - ",celltype))
  lines(density(output_tmpA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="black")
  lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="red")
  lines(density(outputA[nonexpA,]$log2FoldChange,na.rm=TRUE),xlim=c(-3,3),col="salmon")
  text(1.5,1.4,adj=c(0,0),labels="All - baseline",col="black")
  text(1.5,1.35,adj=c(0,0),labels="All - corrected",col="red")
  text(1.5,1.3,adj=c(0,0),labels="NegCtrl - baseline",col="grey")
  text(1.5,1.25,adj=c(0,0),labels="NegCtrl - corrected",col="salmon")
  abline(v=0)
  dev.off()
}


counts <- counts(dds_results, normalized=TRUE)

#Expand IDs that denote duplicate oligos in count/DESeq results
expand_dups<-function(table)
{
  orig_table<-table
  table<-cbind(rownames(table), table)
  colnames(table)[1]<-"Row.names"
  dups<-table[grep("\\(.*\\)$",table$Row.names),]
  dups$Row.names<-gsub("^\\((.*)\\)$","\\1",dups$Row.names)
  final_table<-table[-(grep("\\(.*\\)$",table$Row.names)),] 
  final_table<-final_table[!(is.na(final_table$Row.names)),]
  if(nrow(dups)>0) {
    for(i in 1:nrow(dups))
    {
      dup_id<-unlist(strsplit(dups$Row.names[i],";"))
      dup_exp<-dups[rep(i, length(dup_id)), ]
      dup_exp$Row.names<-dup_id
      final_table<-rbind(dup_exp,final_table)
    }
    rownames(final_table)<-final_table[,1]
    final_table<-final_table[,-1]
  } else {
    final_table<-orig_table
  }
  return(final_table)
}

####TEST code
#ctrl_cols<-c("Plasmid_12","Plasmid_13","Plasmid_14","Plasmid_15","Plasmid_16")
#exp_cols<-c("HTR8_3","HTR8_4","HTR8_5","HTR8_6","HTR8_7")
#Ctrl.Mean=rowMeans(counts[, c(6,7,8,9,10)])
#Exp.Mean=rowMeans(counts[, c(1,2,3,4,5)])
#output_2<-cbind(Ctrl.Mean,Exp.Mean,outputA[,-1])
#output_undup<-expand_dups(output_2)

#tmp_attributeData<-fullattributeData
#counts<-counts
#func_output<-output_undup
#Ctrl.Mean<-Ctrl.Mean
#Exp.Mean<-Exp.Mean
#ctrl_cols<-ctrl_cols
#exp_cols<-exp_cols
####


### Function to perform TTest on individual cell types
CellSpecific_Ttest<-function(tmp_attributeData, counts, func_output, Ctrl.Mean, Exp.Mean, ctrl_cols, exp_cols) {
  
  snp_data<-subset(tmp_attributeData,allele=="ref" | allele=="alt")
  snp_data$comb<-paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,"_",snp_data$haplotype,sep="")
  tmp_ct<-as.data.frame(table(snp_data$comb))
  
  snp_data_pairs<-snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==2,]$Var1,]
  
  snp_data_rejected<-snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq!=2,]$Var1,]
  
  snp_data_ctdata_pairs<-merge(snp_data_pairs,counts,by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  snp_data_ctdata_pairs<-snp_data_ctdata_pairs[order(snp_data_ctdata_pairs$SNP,snp_data_ctdata_pairs$window,snp_data_ctdata_pairs$strand,snp_data_ctdata_pairs$haplotype,snp_data_ctdata_pairs$allele),]
  snp_data_expdata_pairs<-merge(snp_data_pairs,func_output,by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  snp_data_expdata_pairs<-snp_data_expdata_pairs[order(snp_data_expdata_pairs$SNP,snp_data_expdata_pairs$window,snp_data_expdata_pairs$strand,snp_data_expdata_pairs$haplotype,snp_data_expdata_pairs$allele),]
  
  
  ###This is confusing.... but switched the odd/even order here since I am sorting on alt/ref. If ref/alt becomes A/B switch the order here instead of updating the code below. Regardless if this is wrong it will display the A allele in output file correctly.
  evens <- seq(1, nrow(snp_data_pairs), by=2)
  odds <- seq(2, nrow(snp_data_pairs), by=2)
  
  out <- cbind(
    snp_data_expdata_pairs[odds,c(1,2,3,4,5,6,7,9)],
    within(data.frame(
      A.Ctrl.Mean=snp_data_expdata_pairs[odds, "Ctrl.Mean"],
      A.Exp.Mean=snp_data_expdata_pairs[odds, "Exp.Mean"],
      A.log2FC=snp_data_expdata_pairs[odds, "log2FoldChange"],
      A.log2FC_SE=-log10(snp_data_expdata_pairs[odds, "lfcSE"]),
      A.logP=-log10(snp_data_expdata_pairs[odds, "pvalue"]),
      A.logPadj_BH=-log10(snp_data_expdata_pairs[odds, "padj"]),    #BF Correction
      A.logPadj_BF=-log10(snp_data_expdata_pairs[odds, "pvalue"]*(nrow(snp_data_expdata_pairs)/2)),    #BF Correction
      B.Ctrl.Mean=snp_data_expdata_pairs[evens, "Ctrl.Mean"],
      B.Exp.Mean=snp_data_expdata_pairs[evens, "Exp.Mean"],
      B.log2FC=snp_data_expdata_pairs[evens, "log2FoldChange"],
      B.log2FC_SE=-log10(snp_data_expdata_pairs[evens, "lfcSE"]),
      B.logP=-log10(snp_data_expdata_pairs[evens, "pvalue"]),
      B.logPadj_BH=-log10(snp_data_expdata_pairs[evens, "padj"]),    #BF Correction
      B.logPadj_BF=-log10(snp_data_expdata_pairs[evens, "pvalue"]*(nrow(snp_data_expdata_pairs)/2))),{   #BF Correction
        A.logP[is.na(A.logP)] <- 0
        A.logP[A.logP == Inf] <- max(A.logP[is.finite(A.logP)])
        A.logPadj_BH[A.logPadj_BH < 0]<-0
        A.logPadj_BH[A.logPadj_BH == Inf] <- max(A.logPadj_BH[is.finite(A.logPadj_BH)])
        A.logPadj_BF[A.logPadj_BF < 0]<-0
        A.logPadj_BF[A.logPadj_BF == Inf] <- max(A.logPadj_BF[is.finite(A.logPadj_BF)])
        B.logP[is.na(B.logP)] <- 0
        B.logP[B.logP == Inf] <- max(B.logP[is.finite(B.logP)])
        B.logPadj_BH[B.logPadj_BH < 0]<-0
        B.logPadj_BH[B.logPadj_BH == Inf] <- max(B.logPadj_BH[is.finite(B.logPadj_BH)])
        B.logPadj_BF[B.logPadj_BF < 0]<-0
        B.logPadj_BF[B.logPadj_BF == Inf] <- max(B.logPadj_BF[is.finite(B.logPadj_BF)])
      }))
  
  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 |
                        rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) |
                        is.na(rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols])) | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )
  
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1
  
  # t test
  ratios.A <- log((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))
  ratios.B <- log((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  t.pvalue <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(as.numeric(ratios.A[i,]), as.numeric(ratios.B[i,]), var.equal=FALSE, paired=TRUE)$p.value})
  t.stat <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(as.numeric(ratios.A[i,]), as.numeric(ratios.B[i,]), var.equal=FALSE, paired=TRUE)$statistic})
  
  out$LogSkew <- out$B.log2FC - out$A.log2FC
  out$Skew.logP <- ifelse(is.na(t.pvalue), 0, -log10(t.pvalue))
  
  OE_threshold <- -log10(.01)
  is_OE <- out$A.logPadj_BF >= OE_threshold | out$B.logPadj_BF >= OE_threshold
  out$Skew.logFDR <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.logFDR[q_idx] <- -log10(p.adjust(t.pvalue[q_idx],method="BH"))
  
  
  return(out)
}

###############
#### Loop through cell types


full_output<-list()
full_output_var<-list()

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA") next
  message(celltype)
  
  outputA = results(dds_results, contrast=c("condition",celltype,"DNA"))
  
  ctrl_cols<-row.names(colData[colData$condition=="DNA",])
  exp_cols<-row.names(colData[colData$condition==celltype,])
  
  Ctrl.Mean=rowMeans(counts[, colnames(counts) %in% ctrl_cols])
  Exp.Mean=rowMeans(counts[, colnames(counts) %in% exp_cols])
  output_2<-cbind(Ctrl.Mean,Exp.Mean,outputA[,-1])
  output_undup<-expand_dups(output_2)
  
  full_outputA<-merge(fullattributeData, as.matrix(output_undup),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
  full_output[[celltype]]<-full_outputA
  write.table(full_outputA,paste0("results/",file_prefix,"_",celltype,"_",file_date,".out"), row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
  
  outA<-CellSpecific_Ttest(fullattributeData, counts, output_undup, Ctrl.Mean, Exp.Mean, ctrl_cols, exp_cols) 
  full_output_var[[celltype]]<-outA
  write.table(outA,paste0("results/",file_prefix,"_",celltype,"_emVAR_",file_date,".out"), row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}



#### End cell type analysis
###############
###########
##########################
##########################
#########
###### General QC Plots 
#########
counts <- counts(dds_results, normalized=FALSE)
write.table(counts,paste0("results/",file_prefix,"_Counts_",file_date,".out"), row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)
counts <- counts(dds_results, normalized=TRUE)
write.table(counts,paste0("results/",file_prefix,"_NormCounts_",file_date,".out"), row.names=TRUE,col.names=NA,sep="\t",quote=FALSE)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.lm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) 
{   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~ x[ok]), 
           col = col.lm, ...)
}

panel.nlm <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) 
{   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  points(x, y, pch = pch,ylim = ylim, xlim= xlim,col=rgb(144,144,144,75,maxColorValue=255),...)
}

png(file="plots/Cor_mat_log.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,log="xy",pch=16)
dev.off()
png(file="plots/Cor_mat.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,pch=16)
dev.off()
png(file="plots/Cor_mat_2.png",width=3000,height=3000) 
pairs(counts,upper.panel=panel.cor,lower.panel=panel.lm,xlim=c(0,2000),ylim=c(0,2000),pch=16)
dev.off()


ggplot_mpra_scatter<-function(counts, sampleX, sampleY,xmax,ymax) {
  count_df<-as.data.frame(counts)
  ggplot_output<-ggplot(count_df, aes_string(sampleX,sampleY)) +
    theme_minimal() +
    geom_point(alpha = .3,size=1) +
    labs(x=sampleX,y=sampleY) +
    coord_fixed(ratio = 1,xlim=c(0,xmax), ylim=c(0,ymax)) +
    geom_abline(intercept = 0, slope = 1,linetype = 2, size=.75, color=rgb(255,140,0,150,maxColorValue=255))
  return(ggplot_output)
}

xmax<-quantile(counts,.99)
ymax<-quantile(counts,.99)

sampleX<-"Plasmid_r1"
sampleY<-"Plasmid_r2"
ggplot_output<-ggplot_mpra_scatter(counts, sampleX, sampleY,xmax,ymax)
ggsave("plots/plasmid_cor.png",ggplot_output,units="in",width=4,height=4,device="png")

sampleX<-"Plasmid_r1"
sampleY<-"cDNA_r1"
ggplot_output<-ggplot_mpra_scatter(counts, sampleX, sampleY,xmax,ymax)
ggsave("plots/plasmid-Jurkat_cor.png",ggplot_output,units="in",width=4,height=4,device="png")

sampleX<-"cDNA_r1"
sampleY<-"cDNA_r2"
ggplot_output<-ggplot_mpra_scatter(counts, sampleX, sampleY,xmax,ymax)
ggsave("plots/Jurkat_cor.png",ggplot_output,units="in",width=4,height=4,device="png")






plot_logFC<-function(full_output, sample) {
  exp.values<-full_output[full_output$Ctrl.Mean > 10 & !is.na(full_output$Ctrl.Mean),]
  exp.values$Exp.Mean[is.na(exp.values$Exp.Mean)]<-1
  exp.values$log2FoldChange[is.na(exp.values$log2FoldChange)]<-0
  exp.values$padj[is.na(exp.values$padj)]<-1
  exp.values$sig<-"Not Significant"
  exp.values$sig[exp.values$padj <= 0.00001]<-"Active"
  levels(exp.values$sig)<-c("Not Significant", "Active")
  exp.values$sig<-factor(exp.values$sig,levels=c("Not Significant", "Active"))
  
  tmp_plotA<-ggplot(exp.values,aes(x=Ctrl.Mean,y=log2FoldChange,color=sig)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=c("Not Significant"=rgb(0,0,0,200,maxColorValue=255),"Active"=rgb(55,126,184,255,maxColorValue=255))) +
    geom_point(alpha = .3,size=1) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    labs(x="Normalized Tag Count - Plasmids",y=paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .90), 
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) + 
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(255,140,0,150,maxColorValue=255))
  
  cbPalette <- c("#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  tmp_plotB<-1
  tmp_plotB<-ggplot(exp.values,aes(x=Ctrl.Mean,y=log2FoldChange,color=project)) +
    theme_bw() + theme(panel.grid.major = element_line(size = .25,colour = rgb(0,0,0,75,maxColorValue=255)), panel.grid.minor = element_blank()) +
    scale_colour_manual(values=cbPalette) +
    geom_point(alpha = .2,size=1) +
    geom_point(data = subset(exp.values, project == 'negCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp.values, project == 'expCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    geom_point(data = subset(exp.values, project == 'emVarCtrl'), aes(x=Ctrl.Mean,y=log2FoldChange),alpha = .8,size=2) +
    scale_x_log10() +
    #coord_cartesian(xlim = c(10, 1000),ylim = c(-1.5,7.5)) +
    labs(x="Normalized Tag Count - Plasmids",y=paste0(sample," Expression Fold Change log2(RNA/Plasmid)")) +
    theme(legend.position = c(.15, .80), 
          legend.key = element_blank(),
          legend.background = element_rect(color=rgb(0,0,0,150,maxColorValue=255), fill = "white", size = .5, linetype = "solid")) +
    guides(colour = guide_legend(override.aes = list(size=3,alpha=.7), title=NULL)) + 
    geom_abline(intercept = 0, slope = 0,linetype = 1, size=.75, color=rgb(0,0,0,0,maxColorValue=255))
  
  return(list(tmp_plotA,tmp_plotB))
}

for (celltype in levels(colData$condition)) {
  if(celltype=="DNA" | celltype=="EXCLUDE_CELL" ) next
  message(celltype)
  output_tmp<-full_output[[celltype]]
  plot_list<-plot_logFC(output_tmp, celltype)
  ggsave(paste0("plots/logFC_",celltype,".pdf"),plot_list[[1]],units="in",width=8,height=6,device="pdf")
  ggsave(paste0("plots/logFC_",celltype,"_controls.pdf"),plot_list[[2]],units="in",width=8,height=6,device="pdf")
}