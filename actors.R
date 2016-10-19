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

## Plot of the degrees

proportions = table(actor.estimates$Degree.actors)
degree = as.numeric(names(proportions))

proportions = proportions[degree!=0]
degree = degree[degree!=0]

alpha = actor.estimates$est.output[1,1]
slope = -1-alpha

png(
  "./data/actorfull_hollywood.png",
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 6
)
par(mar=c(5,4,1,1)+0.1)
plot(log(degree),log(proportions),
     axes = FALSE,
     ylim = c(0,14))
axis(side=1); axis(side = 2)
lines(log(degree),14+ lgamma(degree-alpha) - lgamma(degree+1) + log(alpha) - lgamma(1-alpha), lty = 2)
dev.off()

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

proportions.er = table(output.er$deg)
degrees.er = as.numeric(names(proportions.er))

png(
  "actor250_erdosrenyi.png",
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 6
)
par(mar=c(5,4,1,1)+0.1)
plot(log(degrees.er),log(proportions.er),axes = FALSE,
     xlab = "log(degree)", ylab = "log(proportions)",
     ylim = c(0,6))
axis(side=1); axis(side = 2)

slope.er = -1-output.er$alpha

abline(a = 6, b = slope.er, lty = 2)

dev.off()

## Fit SBM model

M = as.matrix.network(network.submodel)

#my.model <- BM_bernoulli("SBM_sym",M)

#my.model$estimate()

best.ICL = which(my.model$ICL == max(my.model$ICL))

output.sbm = sbm.powerlaw(my.model)

proportions.sbm = table(output.sbm$deg)
degrees.sbm = as.numeric(names(proportions.sbm))

png(
  "actor250_sbm.png",
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 6
)
par(mar=c(5,4,1,1)+0.1)
plot(log(degrees.sbm),log(proportions.sbm),axes = FALSE,
     xlab = "log(degree)", ylab = "log(proportions)",
     ylim = c(0,6))
axis(side=1); axis(side = 2)

slope.sbm = -1-output.sbm$alpha

abline(a = 6, b = slope.sbm, lty = 2)

dev.off()
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
num.movies2 = 1000
set.seed(8197)
index.subsample2 = sample(1:nr,
                         size=num.movies2,replace=FALSE)
Mact.subsample2 = Mact[index.subsample2,]

fit.subsample2 = model.fit(Mact.subsample2, c(0.5,1))

fit.subsample2$est.output

## Fit ERGM and ER models

number.rows2 = sum(unlist(lapply(
    rowSums(Mact.subsample2!=0),choose.two)))

counter = 1
actor.edge.list2 = matrix(nrow = number.rows2, ncol = 2)

for (i in 1:num.movies2) {
  current.edges = construct.edges(Mact.subsample2[i,])
  if(!any(is.na(current.edges))) {
    actor.edge.list2[counter:(counter-1+ncol(current.edges)),] = t(current.edges)
    counter = counter+ncol(current.edges)
  }
}

unique.actors2 = unique(as.numeric(actor.edge.list2))
actor.edge.list.final2 = 0*actor.edge.list2

for (i in 1:length(unique.actors2)) {
    actor.edge.list.final2[actor.edge.list2
                           == unique.actors2[i]] = i
}

actor.edge.list.final2 =
    uniquecombs(actor.edge.list.final2)

network.submodel2 <- network(actor.edge.list.final2,
                             matrix.type="edgelist",
                             directed = FALSE)

ergm.submodel2 <- ergm(network.submodel2~edges+
                           triangles,
                       control =
                           control.ergm(MCMLE.maxit = 20,
                                        MCMLE.density.guard=exp(4)))

erdos.renyi.submodel2 <- ergm(network.submodel2~edges, control = control.ergm(MCMLE.maxit = 20))

## Comparison

as.numeric(ergm.submodel2$mle.lik)

as.numeric(erdos.renyi.submodel2$mle.lik)

proportions.holly2 = table(fit.subsample2$Degree.actors)[-1]

-fit.subsample2$llik +
    sum(proportions.holly2*
        log(proportions.holly2/sum(proportions.holly2)))

## Erdos-renyi
output.er2 = ergm.powerlaw(erdos.renyi.submodel2)

proportions.er2 = table(output.er2$deg)
degrees.er2 = as.numeric(names(proportions.er2))

png(
  "actor1000_er.png",
  width     = 3.25,
  height    = 3.25,
  units     = "in",
  res       = 1200,
  pointsize = 6
)
par(mar=c(5,4,1,1)+0.1)
plot(log(degrees.er2),log(proportions.er2),axes = FALSE,
     xlab = "log(degree)", ylab = "log(proportions)",
     ylim = c(0,6))
axis(side=1); axis(side = 2)

slope.er2 = -1-output.er2$alpha

abline(a = 6, b = slope.er2, lty = 2)

dev.off()
## Power law

-slope.er2-1

fit.subsample2$est.output[1,1]



