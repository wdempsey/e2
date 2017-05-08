### Functions

dir.proc <- function(num.edges,alpha,theta) {
  n = 1
  deg.vertices = c(1)
  if(runif(1) < (theta+length(deg.vertices)*alpha)/(n+theta)){
    deg.vertices[length(deg.vertices)+1] = 1
    init.edge = c(1,2)
  } else{
    deg.vertices = c(2)
    init.edge = c(1,1)
  }

  edge.set = init.edge
  for (n in 2:(num.edges)) {
    new.edge = vector(length = 2)
    for(i in 1:2) {
      new.node.prob = (theta+length(deg.vertices)*alpha)/(2*(n-1)+i-1+theta)
      if(runif(1) < new.node.prob){
        deg.vertices[length(deg.vertices)+1] = 1
        new.edge[i] = length(deg.vertices)
      } else{
        temp.node = sample(x=1:length(deg.vertices),size=1,prob= (deg.vertices - alpha)/(2*(n-1)+i-1+theta) / (1-new.node.prob))
        deg.vertices[temp.node] = deg.vertices[temp.node]+1
        new.edge[i] = temp.node
      }
    }
    new.edge = new.edge[order(new.edge)]
    edge.set = rbind(edge.set, new.edge)
  }
  return(list("edge.set" = edge.set, "deg.vertices" = deg.vertices))
}

log.lik <- function(data) {
  llik <- function (alpha) {
    return(-sum(lgamma(data-alpha) - lgamma(data+1) + log(alpha) - lgamma(1-alpha)))
  }
  return(llik)
}

alpha.hat <- function(data) {
  optimal = optim(0.5, fn=log.lik(data), lower = 0.001, upper = 0.999)
  return(optimal$par)
}

moment.condition <- function(alpha, num.edges, num.vertices) {
  moment <- function(theta) {
    log.moment = lgamma(theta+1) - log(alpha) - lgamma(theta+alpha) + alpha * log(2 * num.edges)
    num.vertices = length(deg)
    return(abs(log(num.vertices) - log.moment))
  }
  return(moment)
}

sim.proj <- function(alpha, theta, proj.edges, proj.deg) {
  n = nrow(proj.edges)
  deg.vertices = proj.deg
  edge.set = proj.edges[1:nrow(proj.edges), 1:ncol(proj.edges)]
  finish = FALSE
  while (finish != TRUE) {
    n = n + 1
    new.edge = vector(length = 2)
    for(i in 1:2) {
      new.node.prob = (theta+length(deg.vertices)*alpha)/(2*(n-1)+i-1+theta)
      if(runif(1) < new.node.prob){
        finish = TRUE
      } else{
        temp.node = sample(x=1:length(deg.vertices),size=1,prob= (deg.vertices - alpha)/(2*(n-1)+i-1+theta) / (1-new.node.prob))
        deg.vertices[temp.node] = deg.vertices[temp.node]+1
        new.edge[i] = temp.node
      }
    }
    new.edge = new.edge[order(new.edge)]
    if( sum((new.edge[1]==edge.set[,1])*(new.edge[2]==edge.set[,2])) == 0) {
      finish = TRUE
    }
    edge.set = rbind(edge.set, new.edge)
  }
  return(n)
}

avg.degree <- function(num.iters, alpha, theta, proj.edges, proj.deg) {
  scale = nrow(proj.edges)
  total = 0
  for(i in 1:num.iters) {
    total = total + sim.proj(alpha, theta, proj.edges, proj.deg)/scale
  }
  return(total/num.iters * scale)
}

theta.hat <- function(alpha, proj.edge.set, proj.deg) {
  num.iters = 50; max.steps = 25
  num.edges = nrow(proj.edge.set)
  theta.old = optim(1,moment.condition(alpha,num.edges,proj.deg),lower  = -alpha + 0.001, upper = 100)$par
  for(i in 1:max.steps) {
    num.edges = avg.degree(num.iters, alpha, theta.old, proj.edge.set, proj.deg)
    theta.new = optim(1,moment.condition(alpha,num.edges,proj.deg),lower  = -alpha + 0.001, upper = 100)$par
    if( abs((theta.new - theta.old)/theta.old) < 10^(-3)) {break} else{theta.old = theta.new}
  }
  return(theta.new)
}

mg.log.lik <- function(edge.set, deg) {
  num.edges <- nrow(edge.set)
  num.vertices <- length(deg)
  llik <- function(params) {
    alpha = params[1]; theta = params[2];
    return(-1 *(
      (lgamma(num.vertices+theta/alpha) - lgamma(theta/alpha)) -
        (lgamma(2*num.edges+theta) - lgamma(theta)) +
        num.vertices*log(alpha) +
        sum(lgamma(deg[deg>1]-alpha) - lgamma(1-alpha))
    ))
  }
  return(llik)
}

mg.hat <- function(init, edge.set, deg) {
  return(optim(init,mg.log.lik(edge.set, deg), lower = c(0.001,0), upper = c(0.999, 10000) )$par)
}

ks.stat <- function(proportions, alpha) {
  proportions
  degrees = 1:length(proportions)

  log.px = lgamma(degrees-alpha) - lgamma(degrees+1) + log(alpha) - lgamma(1-alpha)

  cdf.px = cumsum(exp(log.px))
  emp.cdf = cumsum(proportions)
  return(max(abs(emp.cdf - cdf.px)/sqrt(cdf.px*(1-cdf.px))))
}

finp.log.lik <- function(V.E, N.sum, deg) {
  llik <- function(params) {
    alpha = params[1]; theta = params[2];
    return(-1 *(
      (lgamma(V.E+theta/alpha) - lgamma(theta/alpha)) -
        (lgamma(N.sum+theta) - lgamma(theta)) +
        V.E*log(alpha) +
        sum(lgamma(deg-alpha) - lgamma(1-alpha))
    ))
  }
  return(llik)
}

finp.hat <- function(init, V.E, N.sum, deg) {
    output = optim(init,finp.log.lik(V.E, N.sum,deg), lower = c(0.001,0), upper = c(0.999, Inf) )
  return(list("params"=output$par,"llik"=output$value))
}

model.fit <- function(data,init) {
  actor.ids = unique(as.vector(data))
  actor.ids = actor.ids[-1]
  actor.ids = actor.ids[order(actor.ids)]

  num.movies.per.actor = rep(0,max(actor.ids))

  for(i in 1:nrow(data)) {
    num.movies.per.actor[data[i,data[i,]!=0]] = num.movies.per.actor[data[i,data[i,]!=0]] +1
  }

  actor.movie.degrees = num.movies.per.actor[num.movies.per.actor!=0]
  num.movies = nrow(data)
  num.roles = sum(actor.movie.degrees)
  num.actors = length(actor.movie.degrees)
  #num.actors = length(unique(as.matrix(Mact))) - 1

  V.E =  num.actors # num.movies
  N.sum = num.roles
  deg = actor.movie.degrees

  finp.output = finp.hat(init, V.E, N.sum, deg)
  finp.params = finp.output$params
  # Calc std err
  library(numDeriv)
  finp.info = hessian(func = finp.log.lik(V.E,N.sum,deg), x = finp.params)
  finp.std.err = sqrt(diag(solve(finp.info)))

  return(list("est.output"=cbind(finp.params, finp.std.err), "N.roles" = num.roles, "N.actors" = num.actors,"Degree.actors"= num.movies.per.actor,"llik"=finp.output$llik))
}

ergm.fit <- function(network.holly) {
    ergm.holly <- ergm(network.holly~edges, control = control.ergm(MCMLE.maxit = 20))
    degree = rowSums(as.matrix(simulate.ergm(ergm.holly)))
    return(list("deg" = degree,
                "alpha" = alpha.hat(degree)))
}

ergm.powerlaw <- function(ergm.holly) {
    degree = rowSums(as.matrix(simulate.ergm(ergm.holly)))
    return(list("deg" = degree,
                "alpha" = alpha.hat(degree)))
}

sbm.fit <- function(network.holly) {
    M = as.matrix.network(network.holly)

    my.model <- BM_bernoulli("SBM_sym",M)

    my.model$estimate()

    best.ICL = which(my.model$ICL == max(my.model$ICL))

    pi = my.model$model_parameters[[best.ICL]]$pi
    membership = my.model$memberships[[best.ICL]]$Z

    M.sim = M*0

    for(i in 1:(nrow(M)-1)) {
        for (j in (i+1):ncol(M)) {
            m.i = sample(1:best.ICL,size=1,prob=membership[i,])
            m.j = sample(1:best.ICL,size=1,prob=membership[j,])
            M.sim[i,j] = rbinom(n = 1, size = 1, prob = pi[m.i, m.j])
        }
    }

    degree = rowSums(M.sim)

    return(alpha.hat(degree))

}

sbm.powerlaw <- function(sbm.model) {

    best.ICL = which(sbm.model$ICL == max(sbm.model$ICL))

    pi = sbm.model$model_parameters[[best.ICL]]$pi
    membership = sbm.model$memberships[[best.ICL]]$Z

    nr = nc = dim(sbm.model$memberships[[best.ICL]]$Z)

    M.sim = matrix(0,nrow = nr, ncol = nc)

    for(i in 1:(nrow(M.sim)-1)) {
        for (j in (i+1):ncol(M.sim)) {
            m.i = sample(1:best.ICL,size=1,prob=membership[i,])
            m.j = sample(1:best.ICL,size=1,prob=membership[j,])
            M.sim[i,j] = rbinom(n = 1, size = 1, prob = pi[m.i, m.j])
        }
    }

    degree = rowSums(M.sim)

    return(list("deg" = degree,
                "alpha" = alpha.hat(degree)))

}


construct.edges <- function(interaction) {
  true.interaction = interaction[interaction!=0]
  true.interaction = true.interaction[order(true.interaction)]
  if(length(true.interaction) >= 2) {
    return(combn(true.interaction,2))
  } else { return(NA)}
}

choose.two <- function(x) {choose(x,2)}


### Finite Case Functions
mg.finite.log.lik <- function(edge.set, deg) {
  num.edges <- nrow(edge.set)
  num.vertices <- length(deg)
  max.vertices = num.vertices
  llik <- function(alpha) {
    theta = -alpha*max.vertices;
    return(-1 *(
      (lgamma(num.vertices+1)) -
        (lgamma(2*num.edges+theta) - lgamma(theta)) +
        num.vertices*log(-alpha) +
        sum(lgamma(deg[deg>1]-alpha) - lgamma(1-alpha))
    ))
  }
  return(llik)
}

mg.finite.hat <- function(init, edge.set, deg) {
  return(optim(init,mg.finite.log.lik(edge.set, deg), lower = -Inf, upper = 0 )$par)
}
