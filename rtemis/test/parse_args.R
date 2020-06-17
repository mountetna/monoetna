test.parse_matrix = function() {
  # a normal matrix with col_names and row_names
  json = '{"matrix":{"col_names":["a", "b", "c"],"row_names":["x","y","z"],"rows":[
    [ 1, 2, 3], [4, 5, 6], [7, 8, 9 ]
  ]}}'

  mat = matrix_arg(fromJSON(json)$matrix)

  checkEquals(colnames(mat), c("a", "b", "c"))
  checkEquals(rownames(mat), c("x", "y", "z"))

  # a matrix with duplicate col_names and row_names
  json = '{"matrix":{"col_names":["a", "a", "c"],"row_names":[null,null,"z"],"rows":[
    [ 1, 2, 3], [4, 5, 6], [7, 8, 9 ]
  ]}}'

  mat = matrix_arg(fromJSON(json)$matrix)

  checkEquals(colnames(mat), c("C1", "C2", "C3"))
  checkEquals(rownames(mat), c("R1", "R2", "R3"))
}
