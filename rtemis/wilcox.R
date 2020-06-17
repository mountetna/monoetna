library(exactRankTests)
library(qvalue)

func_wilcox = function(m1, m2) {
  test_results = lapply(
    1:nrow(m1),
    function(i) {
      test_result = wilcox.exact(
        as.numeric(m1[i,]),
        as.numeric(m2[i,])
      )
      return(
        cbind(
          pvalue=test_result$p.value,
          statistic=test_result$statistic
        )
      )
    }
  )
  results = as.data.frame(t(sapply(test_results,rbind)))
  colnames(results) = c("pvalue", "statistic")
  rownames(results) = rownames(m1)
  tryCatch({
    results$qvalue = qvalue(results$pvalue)$qvalue
  }, error=function(e) {
    # do nothing if we cannot generate a q-value
  })
  return(results)
}
