
func_density = function(values, bw='nrd0') {
  d = density(values, bw)
  return(list( x=d$x, y=d$y ))
}
