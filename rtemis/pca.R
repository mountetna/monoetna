func_pca = function(mat) {
  comps = prcomp(mat)
  return(list(
    rotation=comps$rotation,
    components=comps$x
  ))
}
