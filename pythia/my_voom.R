 library(edgeR)
 library(limma)
 library(methods)
 library(statmod)
 
 run_voom=function( counts, group, num ) {
   y = DGEList(counts=counts,group=group)
   y = y[ rowSums(y$counts) > num, keep.lib.size=FALSE ]
   y <- calcNormFactors(y)
 
   design <- model.matrix(~group)
   v <- voom(y,design)
   fit <- lmFit(v,design)
   fit <- eBayes(fit, robust=TRUE)
   return(fit)
 }
 
 topGenes = function(fit,coef,pval=0.05) {
   return(topTable(fit, coef=coef, number=Inf, sort="p", p=pval))
 }
 
 readCounts = function(counts_file) {
   return(read.delim(counts_file,row.names="gene_id"))
 }
 #args <- commandArgs(trailingOnly = TRUE)
 #counts_table <- args[1]
 #counts = readCounts(counts_table)
 #group = group <- factor(c(2,2,1,2,1,1))
 #fit = run_voom(counts,group)
 #write.table(topGenes(fit,2,1),paste('voom',counts_table,sep='_'),sep='\t')