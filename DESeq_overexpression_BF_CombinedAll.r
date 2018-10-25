rm(list=ls(all=TRUE))
library(DESeq2)
library(svMisc)
library(qvalue)

# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# INPUT=args[1]
# ATTRIBUTES=args[2]
# OUTPUT_PREFIX=args[3]

INPUT <- "/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/Geuv_90K_sep_raw.counts.final"
ATTRIBUTES <- "/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/master.geuvadis.probes.AsiSIremoved.out.class"
OUT_DIR<-"/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/"


analyze <- function(dds, attributes, ctrl_name, ctrl_cols, exp_name, exp_cols) {
  # Subset the data
 
  all_cols <- c(ctrl_cols, exp_cols)
  ctrl_rep <- length(ctrl_cols)
  exp_rep <- length(exp_cols)
  ctrl_range <- 1:ctrl_rep
  exp_range <- (ctrl_rep + 1):(ctrl_rep + exp_rep)
  ncol <- ctrl_rep + exp_rep
  nrow <- nrow(raw.all)
  
  # Filter out zeros
#   uniq_IDs <- apply(attributes, 1, function(x) paste0(x["SNP"], x["Strand"]))
#   # Only delete those where *all* the plasmid counts are 0.
#   zero_idx <- which(apply(raw_counts, 1, function(x) all(x[ctrl_range] == 0) | all(x[exp_range] == 0)))
#   zero_IDs <- uniq_IDs[zero_idx]
#   zero_idx <- which(uniq_IDs %in% zero_IDs)
#   raw_counts <- raw_counts[-zero_idx,]
#   attributes <- attributes[-zero_idx,]
#   message(length(zero_idx), " rows removed (containing zero values)")
#   message("Processing ", dim(raw_counts)[1], " rows.")
  
  # DESeq pipeline
  results <- results(dds, independentFiltering=FALSE, pAdjustMethod="BH",
                    contrast=list(paste0("condition", exp_name), paste0("condition", ctrl_name)))
  
  # Fix the NA pvalues in the results
  results$pvalue[is.na(results$pvalue)] <- 1
  results$padj[is.na(results$padj)] <- 1
  
  # Analysis
  counts <- counts(dds, normalized=TRUE)
  
  odds <- seq(1, nrow(counts), by=2)
  evens <- seq(2, nrow(counts), by=2)
  out <- cbind(
    attributes[odds,],
    within(data.frame(
    A.Ctrl.Mean=rowMeans(counts[odds, ctrl_cols]),
    A.Exp.Mean=rowMeans(counts[odds, exp_cols]),
    A.log2FC=results[odds, "log2FoldChange"],
    A.logP=-log10(results[odds, "pvalue"]),
    A.logPadj=-log10(results[odds, "pvalue"]*(nrow(results)/2)),    #BF Correction
    #   A.logPadj=-log10(results[odds, "padj"]),						#BH Correction
    #   A.logQ=-log10(qvalue(results[odds, "pvalue"])$qvalues),         #q-value correction
    B.Ctrl.Mean=rowMeans(counts[evens, ctrl_cols]),
    B.Exp.Mean=rowMeans(counts[evens, exp_cols]),
    B.log2FC=results[evens, "log2FoldChange"],
    B.logP=-log10(results[evens, "pvalue"]),
    B.logPadj=-log10(results[evens, "pvalue"]*(nrow(results)/2))),{   #BF Correction
    #   B.logPadj=-log10(results[evens, "padj"])),{						#BH Correction
    #   B.logQ=-log10(qvalue(results[evens, "pvalue"])$qvalues)), {		#q-value correction
     A.logP[is.na(A.logP)] <- 0
     A.logP[A.logP == Inf] <- max(A.logP[is.finite(A.logP)])
     A.logPadj[A.logPadj < 0]<-0
     A.logPadj[A.logPadj == Inf] <- max(A.logPadj[is.finite(A.logPadj)])
     #    A.logQ[A.logQ == Inf] <- max(A.logQ[is.finite(A.logQ)])
     B.logP[is.na(B.logP)] <- 0
     B.logP[B.logP == Inf] <- max(B.logP[is.finite(B.logP)])
     B.logPadj[B.logPadj < 0]<-0
     B.logPadj[B.logPadj == Inf] <- max(B.logPadj[is.finite(B.logPadj)])
     #    B.logQ[B.logQ == Inf] <- max(B.logQ[is.finite(B.logQ)])
    }))

  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(counts[odds,ctrl_cols]) < 10 | rowMeans(counts[odds, exp_cols]) < 10 |
                   rowMeans(counts[evens, ctrl_cols]) < 10 | rowMeans(counts[evens, exp_cols]) < 10)
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- counts
  counts1[counts1 == 0] <- 1

  # t test
  ratios.A <- log((counts1[odds, exp_cols]) / rowMeans(counts[odds, ctrl_cols]))
  ratios.B <- log((counts1[evens, exp_cols]) / rowMeans(counts[evens, ctrl_cols]))
  t.pvalue <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(ratios.A[i,], ratios.B[i,], var.equal=FALSE, paired=TRUE)$p.value})
  t.stat <- sapply(1:nrow(ratios.A), function(i) 
    if(i %in% ignore_idx) { NA } else{ t.test(ratios.A[i,], ratios.B[i,], var.equal=FALSE, paired=TRUE)$statistic})
 
  out$LogSkew <- out$B.log2FC - out$A.log2FC
  out$Skew.logP <- ifelse(is.na(t.pvalue), 0, -log10(t.pvalue))
 
 
  message("Checking null T distribution for skew test")
  NOE_threshold <- -log10(.01)
  is_NOE <- out$A.logP < NOE_threshold & out$B.logP < NOE_threshold
  message(length(is_NOE[is_NOE=="TRUE"]))
  NOE_idx<-which(is_NOE)
  pdf(paste(OUT_DIR,"Tnull.",ctrl_name,".",exp_name,".p0.01.pdf",sep=""))
  
  plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10),ylim=c(0,.5),main="T stat distribution")
  x <- seq(-20, 20, length=10000)
  lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255)) 
  dev.off()
  
  NOE_threshold <- -log10(.1)
  pdf(paste(OUT_DIR,"Tnull.",ctrl_name,".",exp_name,".p0.1.pdf",sep=""))
  is_NOE <- out$A.logP < NOE_threshold & out$B.logP < NOE_threshold
  message(length(is_NOE[is_NOE=="TRUE"]))
  NOE_idx<-which(is_NOE)
  plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10),ylim=c(0,.5),main="T stat distribution")
  x <- seq(-20, 20, length=10000)
  lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255)) 
  dev.off()
 
  NOE_threshold <- -log10(.2)
  pdf(paste(OUT_DIR,"Tnull.",ctrl_name,".",exp_name,".p0.2.pdf",sep=""))
  is_NOE <- out$A.logP < NOE_threshold & out$B.logP < NOE_threshold
  NOE_idx<-which(is_NOE)
  message(length(is_NOE[is_NOE=="TRUE"]))
  plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10),ylim=c(0,.5),main="T stat distribution")
  x <- seq(-20, 20, length=10000)
  lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255)) 
  dev.off() 
  
  NOE_threshold <- -log10(.3)
  pdf(paste(OUT_DIR,"Tnull.",ctrl_name,".",exp_name,".p0.3.pdf",sep=""))
  is_NOE <- out$A.logP < NOE_threshold & out$B.logP < NOE_threshold
  NOE_idx<-which(is_NOE)
  message(length(is_NOE[is_NOE=="TRUE"]))
  plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10),ylim=c(0,.5),main="T stat distribution")
  x <- seq(-20, 20, length=10000)
  lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255)) 
  dev.off() 
  
  NOE_threshold <- -log10(.5)
  pdf(paste(OUT_DIR,"Tnull.",ctrl_name,".",exp_name,".p0.5.pdf",sep=""))
  is_NOE <- out$A.logP < NOE_threshold & out$B.logP < NOE_threshold
  NOE_idx<-which(is_NOE)
  message(length(is_NOE[is_NOE=="TRUE"]))
  plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10),ylim=c(0,.5),main="T stat distribution")
  x <- seq(-20, 20, length=10000)
  lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255)) 
  dev.off()
  
# For HepG2 find t-statistic inflation and refit 
#if(exp_name == "HepG2"){
#  message("Refitting t-statistic distribution")
#  NOE_threshold <- -log10(.2)
#  is_NOE <- out$A.logPadj < NOE_threshold & out$B.logPadj < NOE_threshold
#  length(is_NOE[is_NOE=="TRUE"])
#  NOE_idx<-which(is_NOE)
#  #plot(density(t.stat[NOE_idx],na.rm=T),xlim=c(-10,10))
#  #x <- seq(-20, 20, length=10000)
#  #lines(x,dt(x,4), lwd=2, col=rgb(255,0,0,255,maxColorValue=255))
#  z = rank(t.stat[NOE_idx],ties.method = "average");
#  t_z<-qt(z/(length(z)+1),df=(length(exp_cols)-1))
#  #plot(t_z,t.stat[NOE_idx])
#  #abline(0,1)
#  offsets<-t.stat[NOE_idx]/t_z
#  t.stat.cor.factor<-median(offsets,na.rm=T)
#  #hist((offsets))
#  t.stat.cor<-t.stat/t.stat.cor.factor
#  
#  t.pvalue.uncorrected<-t.pvalue
#  t.pvalue<-pt(-abs(t.stat.cor),df=(length(exp_cols)-1))*2
#
#  }
  
  # q-values
  OE_threshold <- -log10(.01)
#  is_OE <- out$A.logQ >= OE_threshold | out$B.logQ >= OE_threshold
  is_OE <- out$A.logPadj >= OE_threshold | out$B.logPadj >= OE_threshold
  out$Skew.logQ <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.logQ[q_idx] <- -log10(qvalue(t.pvalue[q_idx])$qvalues)
	
  OE_threshold <- 0
#  is_OE <- out$A.logQ >= OE_threshold | out$B.logQ >= OE_threshold
  is_OE <- out$A.logPadj >= OE_threshold | out$B.logPadj >= OE_threshold
  out$Skew.logQ.all <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.logQ.all[q_idx] <- -log10(qvalue(t.pvalue[q_idx])$qvalues)

  # BH corrected p-values  
  OE_threshold <- -log10(.01)
#  is_OE <- out$A.logQ >= OE_threshold | out$B.logQ >= OE_threshold
  is_OE <- out$A.logPadj >= OE_threshold | out$B.logPadj >= OE_threshold
  out$Skew.fdr <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.fdr[q_idx] <- -log10(p.adjust(t.pvalue[q_idx],method="BH"))
  
  OE_threshold <- 0
#  is_OE <- out$A.logQ >= OE_threshold | out$B.logQ >= OE_threshold
  is_OE <- out$A.logPadj >= OE_threshold | out$B.logPadj >= OE_threshold
  out$Skew.fdr.all <- rep(NA, dim(out)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t.pvalue)))
  out$Skew.fdr.all[q_idx] <- -log10(p.adjust(t.pvalue[q_idx],method="BH"))


out
}

run_DESeq <- function(raw) {
  condition <- relevel(as.factor(c(rep("Plasmid", 5), rep("12878", 5), rep("19239", 3), rep("HepG2", 5))), "Plasmid")
  colData <- data.frame(row.names=as.character(1:18), condition=condition)
  dds <- DESeqDataSetFromMatrix(countData=raw, colData=colData, design=~condition)
  dds <- DESeq(dds)
  dds
}

# Read in data
raw.all <- read.table(INPUT, header=FALSE, stringsAsFactors=FALSE)
attributes.all <- read.table(ATTRIBUTES, header=TRUE, stringsAsFactors=FALSE)

# Change 1/2 to string properties
attributes.all <- within(attributes.all, {
  Strand <- ifelse(Strand == "1", "pos", "neg")
  Haplotype <- ifelse(Haplotype == "1", "ref", "alt")
})
# Remove the first column of IDs from the raw counts
raw <- raw.all[,-1]

# Subsetting
N <- nrow(raw)
cols.plasmid <- 1:5
cols.12878 <- 6:10
cols.19239 <- 11:13
cols.HepG2 <- 14:18

# Analyze!
dds <- run_DESeq(raw)

NormCounts<-counts(dds, normalized=TRUE)
new_attributes<-cbind(attributes.all[,1:4],rep(c("A","B"),nrow(attributes.all[,1:4])/2))
colnames(new_attributes)[5]<-"allele"
colnames(NormCounts)<-c("plasmid_r1","plasmid_r2","plasmid_r3","plasmid_r4","plasmid_r5","NA12878_r1","NA12878_r2","NA12878_r3","NA12878_r4","NA12878_r5","NA19239_r1","NA19239_r2","NA19239_r3","HepG2_r1","HepG2_r2","HepG2_r3","HepG2_r4","HepG2_r5")
countsTable<-cbind(new_attributes,NormCounts)
write.table(countsTable,file="20150411_SampleCounts.txt",sep="\t",quote=F,row.names=F)

#qvalue
out.12878.2 <- analyze(dds, attributes.all[,1:4], "Plasmid", cols.plasmid, "12878", cols.12878)
out.19239.2 <- analyze(dds, attributes.all[,1:4], "Plasmid", cols.plasmid, "19239", cols.19239)
out.HepG2.2 <- analyze(dds, attributes.all[,1:4], "Plasmid", cols.plasmid, "HepG2", cols.HepG2)
write.table(out.12878.2, "20150102_Geuv.12878.txt", quote=FALSE, row.names=FALSE)
write.table(out.19239.2, "20150102_Geuv.19239.txt", quote=FALSE, row.names=FALSE)
write.table(out.HepG2.2, "20150102_Geuv.HepG2.txt", quote=FALSE, row.names=FALSE)





#Combine pvalues
odds <- seq(1, nrow(attributes.all), by=2)
out.lcl <- cbind(attributes.all[odds,1:4])
  #A Exp    
  is_dir <- (out.12878.2$A.log2FC >= 0 & out.19239.2$A.log2FC >= 0) | (out.12878.2$A.log2FC <= 0 & out.19239.2$A.log2FC <= 0)
  dir_idx<-which(is_dir)
  out.lcl$C.A.ctrl.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.exp.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.log2FC <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.logPadj <- rep(NA, dim(out.lcl)[1])

  out.lcl$C.A.ctrl.mean<- out.12878.2$A.Ctrl.Mean
  out.lcl$C.A.exp.mean<- ((out.12878.2$A.Exp.Mean * 5) + (out.19239.2$A.Exp.Mean * 3)) / 8
  out.lcl$C.A.log2FC<- ((out.12878.2$A.log2FC * 5) + (out.19239.2$A.log2FC * 3)) / 8
  
  is_12878_dir <- (out.12878.2$A.log2FC >= 0 & out.19239.2$A.log2FC < 0 & out.lcl$C.A.log2FC >= 0) | (out.12878.2$A.log2FC <= 0 & out.19239.2$A.log2FC > 0 & out.lcl$C.A.log2FC <= 0)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- (out.12878.2$A.log2FC < 0 & out.19239.2$A.log2FC >= 0 & out.lcl$C.A.log2FC >= 0) | (out.12878.2$A.log2FC > 0 & out.19239.2$A.log2FC <= 0 & out.lcl$C.A.log2FC <= 0)
  dir_19239_idx<-which(is_19239_dir)

  psum<--2*(log(10^-(out.12878.2$A.logP))+log(10^-(out.19239.2$A.logP)))
  fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  out.lcl$C.A.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.19239.2$A.logP)))
  psum<--2*(log(10^-(out.12878.2$A.logP))+log(p_tail))
  out.lcl$C.A.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.12878.2$A.logP)))
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$A.logP)))
  out.lcl$C.A.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))
  
  out.lcl$C.A.logP[out.lcl$C.A.logP == Inf] <- max(out.lcl$C.A.logP[is.finite(out.lcl$C.A.logP)])
  #out.lcl$C.A.logPadj <- -log10(p.adjust((10^-out.lcl$C.A.logP),method="BH"))   #BH Correction
  out.lcl$C.A.logPadj <- -log10(10^-out.lcl$C.A.logP*length(out.lcl$C.A.logP))   #BF Correction
  out.lcl$C.A.logPadj[out.lcl$C.A.logPadj == Inf] <- max(out.lcl$C.A.logPadj[is.finite(out.lcl$C.A.logPadj)])
  out.lcl$C.A.logPadj[out.lcl$C.A.logPadj < 0]<-0

##Old code, does not account for directionality  
  #out.lcl$C.A.log2FC[dir_idx]<- ((out.12878.2$A.log2FC[dir_idx] * 5) + (out.19239.2$A.log2FC[dir_idx] * 3)) / 8
  #psum<--2*(log(10^-(out.12878.2$A.logP))+log(10^-(out.19239.2$A.logP)))
  #fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  #out.lcl$C.A.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  #out.lcl$C.A.logPadj[dir_idx] <- -log10(p.adjust((10^-out.lcl$C.A.logP[dir_idx]),method="BH"))


  #B Exp  
  is_dir <- (out.12878.2$B.log2FC > 0 & out.19239.2$B.log2FC > 0) | (out.12878.2$B.log2FC < 0 & out.19239.2$B.log2FC < 0)
  dir_idx<-which(is_dir)
  out.lcl$C.B.ctrl.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.exp.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.log2FC <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.logPadj <- rep(NA, dim(out.lcl)[1])

  out.lcl$C.B.ctrl.mean<- out.12878.2$B.Ctrl.Mean
  out.lcl$C.B.exp.mean<- ((out.12878.2$B.Exp.Mean * 5) + (out.19239.2$B.Exp.Mean * 3)) / 8
  out.lcl$C.B.log2FC<- ((out.12878.2$B.log2FC * 5) + (out.19239.2$B.log2FC * 3)) / 8
  
  is_12878_dir <- (out.12878.2$B.log2FC >= 0 & out.19239.2$B.log2FC < 0 & out.lcl$C.B.log2FC >= 0) | (out.12878.2$B.log2FC <= 0 & out.19239.2$B.log2FC > 0 & out.lcl$C.B.log2FC <= 0)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- (out.12878.2$B.log2FC < 0 & out.19239.2$B.log2FC >= 0 & out.lcl$C.B.log2FC >= 0) | (out.12878.2$B.log2FC > 0 & out.19239.2$B.log2FC <= 0 & out.lcl$C.B.log2FC <= 0)
  dir_19239_idx<-which(is_19239_dir)

  psum<--2*(log(10^-(out.12878.2$B.logP))+log(10^-(out.19239.2$B.logP)))
  fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  out.lcl$C.B.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.19239.2$B.logP)))
  psum<--2*(log(10^-(out.12878.2$B.logP))+log(p_tail))
  out.lcl$C.B.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.12878.2$B.logP)))
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$B.logP)))
  out.lcl$C.B.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))
  
  out.lcl$C.B.logP[out.lcl$C.B.logP == Inf] <- max(out.lcl$C.B.logP[is.finite(out.lcl$C.B.logP)])
  #out.lcl$C.B.logPadj <- -log10(p.adjust((10^-out.lcl$C.B.logP),method="BH"))  #BH Correction
  out.lcl$C.B.logPadj <- -log10(10^-out.lcl$C.B.logP*length(out.lcl$C.B.logP))   #BF Correction
  out.lcl$C.B.logPadj[out.lcl$C.B.logPadj == Inf] <- max(out.lcl$C.B.logPadj[is.finite(out.lcl$C.B.logPadj)])
  out.lcl$C.B.logPadj[out.lcl$C.B.logPadj < 0]<-0


##Old code, does not account for directionality   
  #out.lcl$C.B.log2FC[dir_idx]<- ((out.12878.2$B.log2FC[dir_idx] * 5) + (out.19239.2$B.log2FC[dir_idx] * 3)) / 8
  #psum<--2*(log(10^-(out.12878.2$B.logP))+log(10^-(out.19239.2$B.logP)))
  #fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  #out.lcl$C.B.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  #out.lcl$C.B.logPadj[dir_idx] <- -log10(p.adjust((10^-out.lcl$C.B.logP[dir_idx]),method="BH"))

  
  #Skew
  OE_threshold <- -log10(0.01)
  is_OE <- out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold
  OE_idx<-which(is_OE)
  is_dir <- ((out.12878.2$LogSkew >= 0 & out.19239.2$LogSkew >= 0) | (out.12878.2$LogSkew <= 0 & out.19239.2$LogSkew <= 0)) & (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold)
  dir_idx<-which(is_dir)
 
  out.lcl$LogSkew.12878 <- rep(NA, dim(out.lcl)[1])
  out.lcl$LogSkew.12878<-out.12878.2$LogSkew
  out.lcl$LogSkew.19239 <- rep(NA, dim(out.lcl)[1])
  out.lcl$LogSkew.19239<-out.19239.2$LogSkew
  out.lcl$LogSkew.Comb<-rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.12878 <- rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.12878<-out.12878.2$Skew.fdr
  out.lcl$Skew.fdr.19239 <- rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.19239<-out.19239.2$Skew.fdr
  
  out.lcl$LogSkew.Comb<- ((out.12878.2$LogSkew * 5) + (out.19239.2$LogSkew * 3)) / 8  
  
  is_12878_dir <- ((out.12878.2$LogSkew >= 0 & out.19239.2$LogSkew < 0 & out.lcl$LogSkew.Comb >= 0) | (out.12878.2$LogSkew <= 0 & out.19239.2$LogSkew > 0 & out.lcl$LogSkew.Comb <= 0)) & (out.lcl$C.A.logPadj >= OE_threshold)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- ((out.12878.2$LogSkew < 0 & out.19239.2$LogSkew >= 0 & out.lcl$LogSkew.Comb >= 0) | (out.12878.2$LogSkew > 0 & out.19239.2$LogSkew <= 0 & out.lcl$LogSkew.Comb <= 0)) & (out.lcl$C.B.logPadj >= OE_threshold)
  dir_19239_idx<-which(is_19239_dir) 
  
  
  out.lcl$C.Skew.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.Skew.fdr <- rep(NA, dim(out.lcl)[1])
  
  psum<--2*(log(10^-(out.12878.2$Skew.logP))+log(10^-(out.19239.2$Skew.logP)))
  out.lcl$C.Skew.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  p_tail<-1-(0.5*(10^-(out.19239.2$Skew.logP)))
  psum<--2*(log(10^-(out.12878.2$Skew.logP))+log(p_tail))
  out.lcl$C.Skew.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  p_tail<-1-(0.5*(10^-(out.12878.2$Skew.logP)))  
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$Skew.logP)))
  out.lcl$C.Skew.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))  

  #calc_fdr_idx<-which(is_dir | is_12878_dir | is_19239_dir)
  calc_fdr_idx<-which(is_dir)

  out.lcl$C.Skew.fdr[calc_fdr_idx] <- -log10(p.adjust((10^-out.lcl$C.Skew.logP[calc_fdr_idx]),method="BH"))

  write.table(out.lcl, "20150102_Geuv.LCL.txt", quote=FALSE, row.names=FALSE)


########################################################################################################################
########################################################################################################################  
######################################################################################################################## 
############### Decided I needed to calculate FDR using all exp sig oligos and not only the dir concordent sites. 

setwd("/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/")
ATTRIBUTES <- "/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/master.geuvadis.probes.AsiSIremoved.out.class"

out.12878.2<-read.table("/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/20150102_Geuv.12878.txt", header=T)
out.19239.2<-read.table("/Users/rtewhey/Dropbox_SP/Functional_Followup/DiffExp_Analysis/analysis_rst_122314/20150102_Geuv.19239.txt", header=T) 
 
attributes.all <- read.table(ATTRIBUTES, header=TRUE, stringsAsFactors=FALSE)

# Change 1/2 to string properties
attributes.all <- within(attributes.all, {
  Strand <- ifelse(Strand == "1", "pos", "neg")
  Haplotype <- ifelse(Haplotype == "1", "ref", "alt")
})
 
  
  #Combine pvalues
odds <- seq(1, nrow(attributes.all), by=2)
out.lcl <- cbind(attributes.all[odds,1:4])
  #A Exp    
  is_dir <- (out.12878.2$A.log2FC >= 0 & out.19239.2$A.log2FC >= 0) | (out.12878.2$A.log2FC <= 0 & out.19239.2$A.log2FC <= 0)
  dir_idx<-which(is_dir)
  out.lcl$C.A.ctrl.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.exp.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.log2FC <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.A.logPadj <- rep(NA, dim(out.lcl)[1])

  out.lcl$C.A.ctrl.mean<- out.12878.2$A.Ctrl.Mean
  out.lcl$C.A.exp.mean<- ((out.12878.2$A.Exp.Mean * 5) + (out.19239.2$A.Exp.Mean * 3)) / 8
  out.lcl$C.A.log2FC<- ((out.12878.2$A.log2FC * 5) + (out.19239.2$A.log2FC * 3)) / 8
  
  is_12878_dir <- (out.12878.2$A.log2FC >= 0 & out.19239.2$A.log2FC < 0 & out.lcl$C.A.log2FC >= 0) | (out.12878.2$A.log2FC <= 0 & out.19239.2$A.log2FC > 0 & out.lcl$C.A.log2FC <= 0)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- (out.12878.2$A.log2FC < 0 & out.19239.2$A.log2FC >= 0 & out.lcl$C.A.log2FC >= 0) | (out.12878.2$A.log2FC > 0 & out.19239.2$A.log2FC <= 0 & out.lcl$C.A.log2FC <= 0)
  dir_19239_idx<-which(is_19239_dir)

  psum<--2*(log(10^-(out.12878.2$A.logP))+log(10^-(out.19239.2$A.logP)))
  fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  out.lcl$C.A.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.19239.2$A.logP)))
  psum<--2*(log(10^-(out.12878.2$A.logP))+log(p_tail))
  out.lcl$C.A.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.12878.2$A.logP)))
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$A.logP)))
  out.lcl$C.A.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))
  
  out.lcl$C.A.logP[out.lcl$C.A.logP == Inf] <- max(out.lcl$C.A.logP[is.finite(out.lcl$C.A.logP)])
  #out.lcl$C.A.logPadj <- -log10(p.adjust((10^-out.lcl$C.A.logP),method="BH"))   #BH Correction
  out.lcl$C.A.logPadj <- -log10(10^-out.lcl$C.A.logP*length(out.lcl$C.A.logP))   #BF Correction
  out.lcl$C.A.logPadj[out.lcl$C.A.logPadj == Inf] <- max(out.lcl$C.A.logPadj[is.finite(out.lcl$C.A.logPadj)])
  out.lcl$C.A.logPadj[out.lcl$C.A.logPadj < 0]<-0

##Old code, does not account for directionality  
  #out.lcl$C.A.log2FC[dir_idx]<- ((out.12878.2$A.log2FC[dir_idx] * 5) + (out.19239.2$A.log2FC[dir_idx] * 3)) / 8
  #psum<--2*(log(10^-(out.12878.2$A.logP))+log(10^-(out.19239.2$A.logP)))
  #fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  #out.lcl$C.A.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  #out.lcl$C.A.logPadj[dir_idx] <- -log10(p.adjust((10^-out.lcl$C.A.logP[dir_idx]),method="BH"))


  #B Exp  
  is_dir <- (out.12878.2$B.log2FC > 0 & out.19239.2$B.log2FC > 0) | (out.12878.2$B.log2FC < 0 & out.19239.2$B.log2FC < 0)
  dir_idx<-which(is_dir)
  out.lcl$C.B.ctrl.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.exp.mean <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.log2FC <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.B.logPadj <- rep(NA, dim(out.lcl)[1])

  out.lcl$C.B.ctrl.mean<- out.12878.2$B.Ctrl.Mean
  out.lcl$C.B.exp.mean<- ((out.12878.2$B.Exp.Mean * 5) + (out.19239.2$B.Exp.Mean * 3)) / 8
  out.lcl$C.B.log2FC<- ((out.12878.2$B.log2FC * 5) + (out.19239.2$B.log2FC * 3)) / 8
  
  is_12878_dir <- (out.12878.2$B.log2FC >= 0 & out.19239.2$B.log2FC < 0 & out.lcl$C.B.log2FC >= 0) | (out.12878.2$B.log2FC <= 0 & out.19239.2$B.log2FC > 0 & out.lcl$C.B.log2FC <= 0)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- (out.12878.2$B.log2FC < 0 & out.19239.2$B.log2FC >= 0 & out.lcl$C.B.log2FC >= 0) | (out.12878.2$B.log2FC > 0 & out.19239.2$B.log2FC <= 0 & out.lcl$C.B.log2FC <= 0)
  dir_19239_idx<-which(is_19239_dir)

  psum<--2*(log(10^-(out.12878.2$B.logP))+log(10^-(out.19239.2$B.logP)))
  fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  out.lcl$C.B.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.19239.2$B.logP)))
  psum<--2*(log(10^-(out.12878.2$B.logP))+log(p_tail))
  out.lcl$C.B.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  
  p_tail<-1-(0.5*(10^-(out.12878.2$B.logP)))
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$B.logP)))
  out.lcl$C.B.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))
  
  out.lcl$C.B.logP[out.lcl$C.B.logP == Inf] <- max(out.lcl$C.B.logP[is.finite(out.lcl$C.B.logP)])
  #out.lcl$C.B.logPadj <- -log10(p.adjust((10^-out.lcl$C.B.logP),method="BH"))  #BH Correction
  out.lcl$C.B.logPadj <- -log10(10^-out.lcl$C.B.logP*length(out.lcl$C.B.logP))   #BF Correction
  out.lcl$C.B.logPadj[out.lcl$C.B.logPadj == Inf] <- max(out.lcl$C.B.logPadj[is.finite(out.lcl$C.B.logPadj)])
  out.lcl$C.B.logPadj[out.lcl$C.B.logPadj < 0]<-0


##Old code, does not account for directionality   
  #out.lcl$C.B.log2FC[dir_idx]<- ((out.12878.2$B.log2FC[dir_idx] * 5) + (out.19239.2$B.log2FC[dir_idx] * 3)) / 8
  #psum<--2*(log(10^-(out.12878.2$B.logP))+log(10^-(out.19239.2$B.logP)))
  #fish_combP<--log10(pchisq(psum,df=4,lower.tail=FALSE))
  #out.lcl$C.B.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  #out.lcl$C.B.logPadj[dir_idx] <- -log10(p.adjust((10^-out.lcl$C.B.logP[dir_idx]),method="BH"))

  
  #Skew
  OE_threshold <- -log10(0.01)
  is_OE <- (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold) & (!is.na(out.lcl$C.A.logPadj) & !is.na(out.lcl$C.B.logPadj))
  OE_idx<-which(is_OE)
  is_dir <- ((out.12878.2$LogSkew >= 0 & out.19239.2$LogSkew >= 0) | (out.12878.2$LogSkew <= 0 & out.19239.2$LogSkew <= 0)) & (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold)
  dir_idx<-which(is_dir)
 
  out.lcl$LogSkew.12878 <- rep(NA, dim(out.lcl)[1])
  out.lcl$LogSkew.12878<-out.12878.2$LogSkew
  out.lcl$LogSkew.19239 <- rep(NA, dim(out.lcl)[1])
  out.lcl$LogSkew.19239<-out.19239.2$LogSkew
  out.lcl$LogSkew.Comb<-rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.12878 <- rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.12878<-out.12878.2$Skew.fdr
  out.lcl$Skew.fdr.19239 <- rep(NA, dim(out.lcl)[1])
  out.lcl$Skew.fdr.19239<-out.19239.2$Skew.fdr
  
  out.lcl$LogSkew.Comb<- ((out.12878.2$LogSkew * 5) + (out.19239.2$LogSkew * 3)) / 8  
  
  is_12878_dir <- ((out.12878.2$LogSkew >= 0 & out.19239.2$LogSkew < 0 & out.lcl$LogSkew.Comb >= 0) | (out.12878.2$LogSkew <= 0 & out.19239.2$LogSkew > 0 & out.lcl$LogSkew.Comb <= 0)) & (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold)
  dir_12878_idx<-which(is_12878_dir)
  is_19239_dir <- ((out.12878.2$LogSkew < 0 & out.19239.2$LogSkew >= 0 & out.lcl$LogSkew.Comb >= 0) | (out.12878.2$LogSkew > 0 & out.19239.2$LogSkew <= 0 & out.lcl$LogSkew.Comb <= 0)) & (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold)
  dir_19239_idx<-which(is_19239_dir) 
  
  
  out.lcl$C.Skew.logP <- rep(NA, dim(out.lcl)[1])
  out.lcl$C.Skew.fdr <- rep(NA, dim(out.lcl)[1])
  
  psum<--2*(log(10^-(out.12878.2$Skew.logP))+log(10^-(out.19239.2$Skew.logP)))
  out.lcl$C.Skew.logP[dir_idx] <- -log10(pchisq(psum[dir_idx],df=4,lower.tail=FALSE))
  p_tail<-1-(0.5*(10^-(out.19239.2$Skew.logP)))
  psum<--2*(log(10^-(out.12878.2$Skew.logP))+log(p_tail))
  out.lcl$C.Skew.logP[dir_12878_idx] <- -log10(pchisq(psum[dir_12878_idx],df=4,lower.tail=FALSE))
  p_tail<-1-(0.5*(10^-(out.12878.2$Skew.logP)))  
  psum<--2*(log(p_tail)+log(10^-(out.19239.2$Skew.logP)))
  out.lcl$C.Skew.logP[dir_19239_idx] <- -log10(pchisq(psum[dir_19239_idx],df=4,lower.tail=FALSE))  

    is_OE_ref <- (out.lcl$C.A.logPadj >= OE_threshold | out.lcl$C.B.logPadj >= OE_threshold) & (!is.na(out.lcl$C.A.logPadj) & !is.na(out.lcl$C.B.logPadj)) & (out.lcl$Haplotype == "ref")

  calc_fdr_idx<-which(is_OE)
  #calc_fdr_idx<-which(is_dir | is_12878_dir | is_19239_dir)
  #calc_fdr_idx<-which(is_dir)

  out.lcl$C.Skew.fdr[calc_fdr_idx] <- -log10(p.adjust((10^-out.lcl$C.Skew.logP[calc_fdr_idx]),method="BH"))

  write.table(out.lcl, file="20150102_Geuv.LCL.txt", quote=FALSE, row.names=FALSE)
  

## Some Plots
# fpm(MPRA_DE,robust=FALSE)
# plot(rowMeans(fpm(MPRA_DE,robust=TRUE)[,1:5]),rowMeans(fpm(MPRA_DE,robust=TRUE)[,6:10]),ylim=c(0,60),xlim=c(0,60),cex=.4,col=rgb(31,31,31,25,maxColorValue=255),xlab="Control",ylab="NA12878",main="MPRA Corrected")
# abline(0,1,col=rgb(255,140,0,200,maxColorValue=255),lty=2)
# 
# plot(rowMeans(counts(MPRA_DE,normalized=TRUE)[,1:5]),rowMeans(counts(MPRA_DE,normalized=TRUE)[,6:10]),cex=.4,col=rgb(31,31,31,25,maxColorValue=255),xlab="Control",ylab="NA12878",main="MPRA Corrected")
# abline(0,1,col=rgb(255,140,0,200,maxColorValue=255),lty=2)

## end plots


# OUTFILE=paste(OUTPUT_PREFIX,".results.txt",sep="")
# write.table(out,file=OUTFILE,quote=F,col.names=NA)
