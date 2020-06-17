func_center = function(mat) {
  # recenter a matrix
  return(t(scale(t(mat), scale=FALSE)))
}

func_scale = function(mat) {
  return(t(scale(t(mat), center=FALSE)))
}

func_log = function(mat) {
  return(log(mat))
}

func_transpose = function(mat) {
  return(
    t(mat)
  )
}

func_correlation = function(m1,m2=NULL,method="spearman") {
  return(
    cor(m1, m2, method=method)
  )
}
