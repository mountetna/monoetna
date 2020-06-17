library(beeswarm)

func_beeswarm = function(data) {
  pdf('/dev/null')
  swarm = beeswarm(data, do.plot=F)
  dev.off()
  return(swarm[ , c("x", "y")])
}
