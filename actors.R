## April 8, 2015
## Actor network data set from
## http://www3.nd.edu/~networks/resources/actor/actor.dat.gz

## Actor network data: (based on www.imdb.com)
##
## Gzipped ASCII file (3.9MB)
##
## Each line corresponds to one movie
## Each number represents actor:
##
## Data-format:
## number_1 number_2 ... number_k    (k actors who play in the same movie)
##
## For more information, please contact hjeong@nd.edu

source("./functions.R")

act <- read.csv("./data/actors.txt",header=FALSE,sep=" ")

nr <- nrow(act)
nc <- ncol(act)

Mact <- matrix(0,nr,nc)

for(i in 1:nr){for(j in 1:nc){
	if(is.na(act[i,j])==TRUE){Mact[i,j] <- -1}else{Mact[i,j] <- act[i,j]}
}}

Mact <- Mact+1

actor.estimates = model.fit(Mact,c(0.5,1))

actor.estimates$est.output
actor.estimates$N.roles
actor.estimates$N.actors
actor.estimates$llik

alpha.hat(actor.estimates$Degree.actors)

## Plot of the degrees

proportions = table(actor.estimates$Degree.actors)
degree = as.numeric(names(proportions))

proportions = proportions[degree!=0]
degree = degree[degree!=0]

alpha = actor.estimates$est.output[1,1]
slope = -1-alpha

png(
  "deg_distn_graphs.png",
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 6
)
par(mar=c(3,3.5,1,1)+0.1, mfrow = c(2,2))
plot(log(degree),log(proportions),
     axes = FALSE,
     xlab = "", ylab = "",
     ylim = c(0,12),
     pch = "*")
axis(side=1, cex.axis = 0.75)
axis(side = 2,cex.axis = 0.75)
mtext("log(degree)", side=1, line = 2,cex =0.75)
mtext("log(proportions)", side=2, line = 2, cex = 0.75)
text(5,11,"(A)")
#dev.off()

## Construct a small subsample of movies
library(mgcv); library(blockmodels)
library(statnet)
num.movies = 250
set.seed(8197)
index.subsample = sample(1:nr,
                         size=num.movies,replace=FALSE)
Mact.subsample = Mact[index.subsample,]

fit.subsample = model.fit(Mact.subsample, c(0.5,1))

fit.subsample$est.output

alpha = fit.subsample$est.output[1,1]

lines(log(degree),12.75+ lgamma(degree-alpha) - lgamma(degree+1) + log(alpha) - lgamma(1-alpha), lty = 2)

## Fit ERGM and SBM models

number.rows = sum(unlist(lapply(rowSums(Mact.subsample!=0),choose.two)))

counter = 1; actor.edge.list = matrix(nrow = number.rows, ncol = 2)

for (i in 1:num.movies) {
  current.edges = construct.edges(Mact.subsample[i,])
  if(!any(is.na(current.edges))) {
    actor.edge.list[counter:(counter-1+ncol(current.edges)),] = t(current.edges)
    counter = counter+ncol(current.edges)
  }
}

unique.actors = unique(as.numeric(actor.edge.list))
actor.edge.list.final = 0*actor.edge.list

for (i in 1:length(unique.actors)) {
  actor.edge.list.final[actor.edge.list == unique.actors[i]] = i
}

actor.edge.list.final = uniquecombs(actor.edge.list.final)

network.submodel <- network(actor.edge.list.final,matrix.type="edgelist", directed = FALSE)

erdos.renyi.submodel <- ergm(network.submodel~edges, control = control.ergm(MCMLE.maxit = 20))

output.er = ergm.fit(network.submodel)

proportions.true = degreedist(network.submodel)
degrees.true = c(1:13,15:18,20)

proportions.er = table(output.er$deg)[-1]
degrees.er = as.numeric(names(proportions.er))

plot(log(degrees.true),log(proportions.true),axes = FALSE,
     xlab = "", ylab = "",
     ylim = c(0,6),
     pch = "*")
axis(side=1, cex.axis = 0.75)
axis(side = 2,cex.axis = 0.75)
mtext("log(degree)", side=1, line = 2,cex =0.75)
mtext("log(proportions)", side=2, line = 2, cex = 0.75)
lines(loess(log(proportions.er)~log(degrees.er)), lty = 2)
text(2.75,5.5,"(B)")

slope.er = -1-output.er$alpha

abline(a = 6, b = slope.er, lty = 2)

## Fit SBM model

M = as.matrix.network(network.submodel)

##my.model <- BM_bernoulli("SBM_sym",M)

##my.model$estimate()

load("sbm_model.RData")

best.ICL = which(my.model$ICL == max(my.model$ICL))

output.sbm = sbm.powerlaw(my.model)

proportions.sbm = table(output.sbm$deg)[-1]
degrees.sbm = as.numeric(names(proportions.sbm))

plot(log(degrees.true),log(proportions.true),axes = FALSE,
     xlab = "", ylab = "",
     ylim = c(0,6),
     pch = "*")
lines(loess(log(proportions.sbm)~log(degrees.sbm)), lty = 2)
axis(side=1, cex.axis = 0.75)
axis(side = 2,cex.axis = 0.75)
mtext("log(degree)", side=1, line = 2,cex =0.75)
mtext("log(proportions)", side=2, line = 2, cex = 0.75)
text(2.75,5.5,"(C)")

slope.sbm = -1-output.sbm$alpha

abline(a = 6, b = slope.sbm, lty = 2)

#dev.off()
### Comparison of all
### three methods.

## Log.likelihood

max(my.model$ICL)

as.numeric(erdos.renyi.submodel$mle.lik)

proportions.holly = table(fit.subsample$Degree.actors)[-1]

-fit.subsample$llik +
    sum(proportions.holly*
        log(proportions.holly/sum(proportions.holly)))

## Power law

-slope.sbm-1

-slope.er-1

fit.subsample$est.output[1,1]

## Construct
## larger subsamples
## of movies.  Only fit
## our model and ERGM
num.samples = 1000; set.seed(5131)
movie.sample = sample(1:nrow(Mact),size = num.samples,replace = FALSE)

test.data = Mact[movie.sample,]

number.rows = sum(unlist(lapply(rowSums(test.data!=0),choose.two)))

counter = 1; actor.edge.list = matrix(nrow = number.rows, ncol = 2)

for (i in 1:num.samples) {
  current.edges = construct.edges(test.data[i,])
  if(!any(is.na(current.edges))) {
    actor.edge.list[counter:(counter-1+ncol(current.edges)),] = t(current.edges)
    counter = counter+ncol(current.edges)
  }
}

unique.actors = unique(as.numeric(actor.edge.list))
actor.edge.list.final = 0*actor.edge.list

for (i in 1:length(unique.actors)) {
  actor.edge.list.final[actor.edge.list == unique.actors[i]] = i
}

actor.edge.list.final = uniquecombs(actor.edge.list.final)

# attr(actor.edge.list,"vnames") = unique.actors
# attr(actor.edge.list,"n") = length(unique.actors)

network.test <- network(actor.edge.list.final,matrix.type="edgelist", directed = FALSE)

ergm.model2 <- ergm(network.test~edges+triangles, control = control.ergm(MCMLE.maxit = 20))

er.model2 <- ergm(network.test~edges, control = control.ergm(MCMLE.maxit = 20))

proportions.test = degreedist(network.test)
degrees.test = c(1:29,32,34,36,49)

deg = degreedist(simulate.ergm(ergm.model2))

optim(0.5,llik(deg),lower = 0.001,upper = 0.999)$par

## Comparison

as.numeric(er.model2$mle.lik)
as.numeric(ergm.model2$mle.lik)

proportions.holly2 = table(fit.subsample2$Degree.actors)[-1]

-fit.subsample2$llik +
    sum(proportions.holly2*
        log(proportions.holly2/sum(proportions.holly2)))

## Erdos-renyi
output.ergm2 = ergm.powerlaw(ergm.model2)

proportions.ergm2 = table(output.ergm2$deg)[-1]
degrees.ergm2 = as.numeric(names(proportions.ergm2))

output.er2 = ergm.powerlaw(er.model2)

proportions.er2 = table(output.er2$deg)[-1]
degrees.er2 = as.numeric(names(proportions.er2))


plot(log(degrees.test),log(proportions.test),axes = FALSE,
     xlab = "", ylab = "",
     ylim = c(0,8),
     pch = "*")
axis(side=1, cex.axis = 0.75)
axis(side = 2,cex.axis = 0.75)
mtext("log(degree)", side=1, line = 2,cex =0.75)
mtext("log(proportions)", side=2, line = 2, cex = 0.75)
lines(loess(log(proportions.er2)~
                log(degrees.er2),
            span = 0.4), lty = 2)
lines(loess(log(proportions.ergm2)~
                log(degrees.ergm2),
            span = 0.4), lty = 2)#, col = "red")
text(3.5,7.5,"(D)")
dev.off()

slope.er2 = -1-output.er2$alpha

abline(a = 6, b = slope.er2, lty = 2)

## Power law

-slope.er2-1

fit.subsample2$est.output[1,1]



