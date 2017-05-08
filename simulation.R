## Take output of E2 estimation
## and compute power law under
## simulation from the Hollywood
## model.
source("./functions.R")

## Simulate 2000 edges from a
## hollywood model
library(mgcv); library(blockmodels);
library(statnet)
set.seed(4092)
num.edges = 2000
alpha = 0.6588141
theta = 4.2139409
num.iters = 10
fit.powerlaw =
    fit.proj.powerlaw = vector(length = num.iters)
fit.theta =
    fit.proj.theta = vector(length = num.iters)
fit.ergm.powerlaw =
    fit.sbm.powerlaw = vector(length = num.iters)

for(i in 1:num.iters) {
    ## E2 and projected E2 fit
    sim.holly = dir.proc(num.edges, alpha, theta)
    proj.edges = uniquecombs(sim.holly$edge.set)
    fit.holly = model.fit(sim.holly$edge.set)
    fit.projholly = model.fit(proj.edges)
    fit.powerlaw[i] = fit.holly$est.output[1,1]
    fit.proj.powerlaw[i] = fit.projholly$est.output[1,1]
    fit.theta[i] = fit.holly$est.output[2,1]
    fit.proj.theta[i] = fit.projholly$est.output[2,1]

    projnetwork.holly <- network(proj.edges,Matrix.type="edgelist", directed = FALSE)
    fit.ergm.powerlaw[i] = ergm.fit(projnetwork.holly)
    fit.sbm.powerlaw[i] = sbm.fit(projnetwork.holly)
}

mean(fit.powerlaw); sd(fit.powerlaw)
mean(fit.proj.powerlaw); sd(fit.proj.powerlaw)
mean(fit.ergm.powerlaw); sd(fit.ergm.powerlaw)
mean(fit.sbm.powerlaw); sd(fit.sbm.powerlaw)

mean(fit.theta); sd(fit.theta)
mean(fit.proj.theta); sd(fit.proj.theta)


## Simulate 2000 edges from a
## hollywood model.
## Check Standard Errors
library(mgcv); library(blockmodels);
library(statnet)
set.seed(4091)
num.edges = 2000
seq.alpha = c(0.2,0.4,0.6)
seq.theta = c(5,10,20)
num.iters = 100
fit.alpha =
    fit.theta = vector(length = num.iters)
fit.sd.alpha =
  fit.sd.theta = vector(length = num.iters)
complete.results = rep(0,0)

for (alpha in seq.alpha) {
  for (theta in seq.theta) {
    for(i in 1:num.iters) {
      ## E2 and projected E2 fit
      sim.holly = dir.proc(num.edges, alpha, theta)
      fit.holly = model.fit(sim.holly$edge.set, c(0.5,1))
      fit.alpha[i] = fit.holly$est.output[1,1]
      fit.theta[i] = fit.holly$est.output[2,1]
      fit.sd.alpha[i] = fit.holly$est.output[1,2]
      fit.sd.theta[i] = fit.holly$est.output[2,2]
    }
    
    sim.output = matrix(
      c(alpha, mean(fit.alpha),sd(fit.alpha), mean(fit.sd.alpha),
        theta, mean(fit.theta), sd(fit.theta), mean(fit.sd.theta)),
      nrow = 2, ncol = 4, byrow = TRUE)
    complete.results = rbind(complete.results,sim.output)
    print(paste("Finished with alpha =",alpha,"and theta =",theta))
  }
}


### Histograms of hat.alpha and hat.theta stacked for alpha =0.6, theta = 10
set.seed(4091)
num.edges = 10000
alpha = 0.6
theta = 5
fit.alpha =
  fit.theta = vector(length = num.iters)
fit.sd.alpha =
  fit.sd.theta = vector(length = num.iters)

for(i in 1:num.iters) {
  ## E2 and projected E2 fit
  sim.holly = dir.proc(num.edges, alpha, theta)
  fit.holly = model.fit(sim.holly$edge.set, c(0.5,1))
  fit.alpha[i] = fit.holly$est.output[1,1]
  fit.theta[i] = fit.holly$est.output[2,1]
  fit.sd.alpha[i] = fit.holly$est.output[1,2]
  fit.sd.theta[i] = fit.holly$est.output[2,2]
}


png("MLE_hist_06_5_10000.png", width = 6.5, height = 3, units = "in", res = 300)
par(mfrow = c(1,2) , mar = c(5,4,1,1)+0.1)
hist(fit.alpha, breaks = 10, xlab = expression(paste("MLE for ", alpha)), ylab = "Frequency", main = "" )
abline(v = alpha, lty = 2, lwd = 3)
hist(fit.theta, breaks = 10, xlab = expression(paste("MLE for ", theta)), ylab = "Frequency", main = "" )
abline(v = theta, lty = 2, lwd = 3)
dev.off()
