 library(edgeR)
 library(limma)
 library(methods)
 library(statmod)
 
 func_diff_exp=function( counts, group, pval ) {
   non_zero_columns = which(colSums(counts) > 0)
   non_zero_rows = which(rowSums(counts) > 0)
   counts = counts[ non_zero_rows, non_zero_columns]
   group = group[ non_zero_columns ]
   num_counts = ncol(counts)*10
   y = DGEList(counts=counts,group=group)
   y = y[ rowSums(y$counts) > num_counts, keep.lib.size=FALSE ]
   y = calcNormFactors(y)
 
   design = model.matrix(~group)
   v = voom(y,design)
   fit = lmFit(v,design)
   fit = eBayes(fit, robust=TRUE)
   return(topTable(fit, coef=2, number=Inf, sort="p", p=pval))
}
