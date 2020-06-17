# input transformations

matrix_arg = function(mat) {
  m = as.data.frame(t(as.data.frame(mat$rows)))
  rownames(m) = if (length(mat$row_names) != length(unique(mat$row_names))) paste('R', 1:length(mat$row_names),sep="") else mat$row_names
  colnames(m) = if (length(mat$col_names) != length(unique(mat$col_names))) paste('C', 1:length(mat$col_names),sep="") else mat$col_names
  return(m)
}

vector_arg = function(vec) {
  rvec = unlist(lapply(vec, function(v) return(if (is.null(v$value)) NA else v$value)))
  names(rvec) = unlist(lapply(vec, function(v) return(if (is.null(v$label)) NA else v$label)))
  return(rvec)
}

# the main input transformation loop
parse_args = function(args) {
  return(
    lapply(args, function(arg) {
      if (class(arg) == "list" && "matrix" %in% names(arg)) return(matrix_arg(arg$matrix))
      if (class(arg) == "list" && "vector" %in% names(arg)) return(vector_arg(arg$vector))
      return(arg)
    })
  )
}

vector_output = function(vec) {
  return(
    list(
      vector=lapply(
        seq_along(vec), function(i) {
          return(
            list(
              label=names(vec)[[i]],
              value=transform_outputs(vec[[i]])
            )
          )
        }
      )
    )
  )
}

# output transformations
matrix_output = function(mat) {
  row_names = if (is.null(rownames(mat))) paste('R', 1:nrow(mat),sep="") else rownames(mat)
  col_names = if (is.null(colnames(mat))) paste('C', 1:ncol(mat),sep="") else colnames(mat)
  return(
    list(
      matrix = list(
        row_names = as.list(row_names),
        col_names = as.list(col_names),
        rows = array(
          lapply(split(as.matrix(mat),1:nrow(mat)), as.list)
        )
      )
    )
  )
}

# the main output transformation loop

transform_outputs = function(output) {
  if (class(output) == "matrix" || class(output) == "data.frame") return(matrix_output(output))
  if (class(output) == "list") return(vector_output(output))

  return(output)
}
