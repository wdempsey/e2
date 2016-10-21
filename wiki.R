#######################################################
source("functions.R")

wiki <- read.csv("./data/wiki.txt",header=FALSE, sep=" ")

nr.wiki <- nrow(wiki)
nc.wiki <- ncol(wiki)

mx.wiki <- max(wiki)
M.wiki <- rep(0,mx.wiki)

for(i in 1:nr.wiki){
	M.wiki[wiki[i,1]] <- M.wiki[wiki[i,1]]+1;
}

zero.wiki <- c()
for(i in 1:length(M.wiki)){
	if(M.wiki[i]==0){zero.wiki <- c(zero.wiki,i)}
}

dmz.wiki <- M.wiki[-zero.wiki]

wiki.params = mg.hat(c(0.5,1), edge.set=wiki,deg=dmz.wiki)
# Calc the std dev
library(numDeriv)
wiki.info = hessian(func=mg.log.lik(edge.set=wiki, deg = dmz.wiki), x = wiki.params)
std.err = sqrt(diag(solve(wiki.info)))

round(cbind(wiki.params, std.err),3)

alpha.hat(data=dmz.wiki)

