source("./functions.R")
library(numDeriv)

airport <- read.csv(url("http://opsahl.co.uk/tnet/datasets/USairport_2010.txt"),header=FALSE,sep=" ")

num.edges = sum(airport[,3])
num.vertices = length(unique(c(airport[,1],airport[,2])))
proj.deg = deg = vector(length = num.vertices)

for(i in 1:num.vertices) {
  deg[i] =  sum(airport[airport[,1]==i | airport[,2] == i,3])
  proj.deg[i] =  length(airport[airport[,1]==i | airport[,2] == i,3])
}

airport.log.lik <- function(num.vertices,num.edges, deg) {
  llik <- function(params) {
      alpha = params[1]
      if(alpha > 0) {
          theta = params[2]
          return(-1 *(
              (lgamma(num.vertices+theta/alpha) - lgamma(theta/alpha)) -
              (lgamma(2*num.edges+theta) - lgamma(theta)) +
              num.vertices*log(alpha) +
              sum(lgamma(deg[deg>1]-alpha) - lgamma(1-alpha))
          ))
      } else {
          theta = - num.vertices * alpha
          return(-1 *(
              num.vertices*log(-alpha) +
              - (lgamma(2*num.edges+theta) - lgamma(theta)) +
              sum(lgamma(deg[deg>1]-alpha) - lgamma(1-alpha))
          ))
      }

  }
  return(llik)
}

airport.hat.positive <- function(init, num.vertices, num.edges, deg) {
  return(optim(init,airport.log.lik(num.vertices,num.edges, deg), lower = c(0.005,0.001), upper = c(0.99, 1000000) ))
}

airport.hat.negative <- function(init, num.vertices, num.edges, deg) {
  return(optim(init,airport.log.lik(num.vertices,num.edges, deg), upper = c(-0.00001) ))
}


## Inf pop case: 0 < \alpha < 1
init= c(0.5,1)
positive.output = airport.hat.positive(init,num.vertices,num.edges,deg)
pos.airport.info = hessian(func=airport.log.lik(num.vertices,num.edges, deg), x=positive.output$par)
pos.airport.stderr = sqrt(diag(solve(pos.airport.info)))
cbind(positive.output$par,pos.airport.stderr)
positive.output$value/10^9

## Comparison to Yule Law
alpha.hat(deg)


## Finite case:
init = -5
negative.output = airport.hat.negative(init,num.vertices,num.edges,deg)
neg.airport.info = hessian(func=airport.log.lik(num.vertices,num.edges, deg), x=negative.output$par)
neg.airport.stderr = sqrt(diag(solve(neg.airport.info)))
cbind(negative.output$par,neg.airport.stderr)
negative.output$value/10^9

-positive.output$value/10^9 > -negative.output$value/10^9
