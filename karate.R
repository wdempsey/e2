karate <- read.csv("./data/zachary-weighted.txt",header=FALSE,sep=" ")
karate <- as.matrix(karate)

source("./functions.R")

### Construct Edge Set
deg = rowSums(karate)
edge.set = matrix(nrow = 1, ncol = 2)

for(i in 1:nrow(karate)) {
  for(j in i:ncol(karate)) {
    edge.set = rbind(edge.set,cbind(rep(i,karate[i,j]),rep(j,karate[i,j]) ))
  }
}
edge.set = edge.set[-1,]

karate.alpha.hat = mg.finite.hat(init=-0.5,edge.set=edge.set, deg=deg)
library(numDeriv)
karate.info = hessian(func = mg.finite.log.lik(edge.set,deg), x = karate.alpha.hat)
karate.std.err = sqrt(diag(solve(karate.info)))

cbind(karate.alpha.hat, karate.std.err)

theta.hat = -karate.alpha.hat*length(deg)

theta.std.err = karate.std.err * length(deg)

cbind(theta.hat, theta.std.err)
