# Libs ----
library(MASS)

# Config ----
if (F){
  setwd("~/Fysmat/8 Semester V2022/Ber. Stat. Met/TMA4300/Exe2bsm")
  options(error=recover)
  options(error=NULL)
}
# Problem 1
# e) Implementation of MCMC

# Functions ----
targetdist = function(sigma_u2, tau, Tval){
  lnp = log(sigma_u2) - (Tval-1)*tau %*% Q %*% tau / (2*sigma_u2)
  # lnp = -(Tval-1)*log(sigma_u2) - (Tval-1)*tau %*% Q %*% tau / (2*sigma_u2)
  return(exp(lnp))
}


link = function(tau){
  return(exp(tau)/(1+exp(tau)))
}

logbinom = function(n, y, tau){
  return(y*tau - n*log(1+exp(tau)))
}

acceptRatio = function(n, y, tauProp, tau){
  return(exp(y*(tauProp - tau) + n*log((1+exp(tau))/(1+exp(tauProp)))))
}


mhRW = function(tau, sigma, yt, t, normVec=NA){
  # separate this maybe? as jostein did
  if (t==1){
    mu_ab = tau[2]
    sigma_aa = sigma
  }
  else if (t==366){
    mu_ab = tau[365]
    sigma_aa = sigma
  }
  else{
    mu_ab = 1/2 * (tau[t-1] + tau[t+1])
    sigma_aa = sigma/2
  }
  # prop_tau = rnorm(1,mean=mu_ab, sd=sqrt(sigma_aa))
  prop_tau = normVec[t]*sigma + mu_ab
  n = ifelse(t==60,10,39)
  # if (t == 60){
    # browser()
  # }
  # ratio = dbinom(yt, n, link(prop_tau)) / dbinom(yt, n, link(tau[t]))
  # ratio = exp(logbinom(n, yt, prop_tau)/logbinom(n, yt, tau[t]))
  ratio = acceptRatio(n, yt, prop_tau, tau[t])
  if (runif(1) < min(c(1,ratio))){
    return(list(tau=prop_tau, accepted=1))
  }
  else{return(list(tau=tau[t], accepted=0))}
}

# mhRW = function(tau, sigma, yt, t){
#   mu = 0 - 1/Q[t,t] * Q[t,-t] %*% tau[-t]
#   prop_tau = rnorm(n=1,mu=mu, Sigma=sigma/Q[t,t])
#   n = ifelse(t==60,10,39)
#   ratio = dbinom(yt, n, link(prop_tau)) / dbinom(yt, n, link(tau[t]))
#   # if (runif(1) < accept(n, yt, prop_tau, tau[t])){
#   if (runif(1) < min(c(1,ratio))){
#     return(list(tau=prop_tau, accepted=1))
#   }
#   else{return(list(tau=tau[t], accepted=1))}
# }

gibbs = function(tau, Q, alpha = 2, beta = 0.05){
  g = rgamma(1, alpha, t(tau) %*% Q %*% tau + beta)
  return(1/g)
}

mcmcRW = function(N, dt, sigma0=0.1){
  # Allocate memory
  # browser()
  tau = matrix(NA, nrow=N, ncol = 366)
  sigma = numeric(length = N)
  tau_i = numeric(length = 366)
  normMat = matrix(rep(rnorm(366), N), nrow = N, ncol=366)
  # bc = choose(dt$n.years, dt$n.rain) # Binom coeff. not needed?!
  
  # Find init vals
  pi0 = dt$n.rain/dt$n.years
  # tau[1,] = log(pi0/(1-pi0)) # empiric tau
  tau[1,] = rnorm(366) # tau Drawn from normal distr.
  sigma[1] = sigma0
  # Make Q matrix
  Q=matrix(0, nrow = 366, ncol = 366)
  diag(Q)=2
  Q[c(1, length(Q))]=1
  Q[abs(row(Q) - col(Q)) == 1] <- -1
  # Run mcmc for N iterations
  accepted = 0
  for (i in 2:N){
    tau_i = tau[i-1,]
    sigma_i = sigma[i-1]
    for (t in 1:366){
      # rtemp= mhRW(tau_i, sigma_i, dt$n.rain[t], t)
      rtemp= mhRW(tau_i, sqrt(sigma_i), dt$n.rain[t], t, normMat[i,])
      tau[i,t] = rtemp$tau
      accepted = accepted + rtemp$accepted
    }
    # sigma[i] = gibbs(tau_i, Q)
    # sigma[i] = 1/rgamma(1,2,0.05) # Non gibbs
    
    # tQtn = t(tau[i,]) %*% Q %*% tau[i,]
    tQtn = sum((tau[i,-366] - tau[i,-1])^2)
    if (is.infinite(tQtn)){
      browser()
      tQtn=0
    }
    sigma[i] = 1/rgamma(1, 2 + (366-1)/2, 0.05 + 0.5*tQtn) # Gibbs inline
    # tQto = t(tau_i) %*% Q %*% tau_i
    # sigma[i] = 1/rgamma(1, 2 + (366-1)/2, 0.05 + 0.5*tQto) # Gibbs inline
    # sigma[i] = sigma0
    if (i%%(N/10)==0){
      print(i/N*100)
      print(accepted/(i*366))
    }
  }
  return(list(tau=tau, sigma=sigma))
}



## setup ----

load("./rain.rda")
# rain$tau_t
rain$pi = rain$n.rain/rain$n.years
rain$tau_t= log((rain$pi)/(1-rain$pi))
rain$tau_t[60]=0

Tval = length(rain$tau_t)
Q=matrix(0, nrow = Tval, ncol = Tval)
diag(Q)=2
Q[c(1, length(Q))]=1

Q[abs(row(Q) - col(Q)) == 1] <- -1
# Q[1:5, 1:5]
# Q[361:366, 361:366]

sigma_u2 = 0.01


tau = rain$tau_t
lnp = log(sigma_u2) - (Tval-1)*tau %*% Q %*% tau / (2*sigma_u2)
target = targetdist(sigma_u2, rain$tau_t, Tval)
# target

# tests ----
# mhRW(rain$tau_t, 0.1, rain$n.rain[1], 1)

set.seed(321)
results = mcmcRW(10000, rain)
sum(is.na(results$tau))
sum(is.na(results$sigma))

par(mfrow=c(3,3))
if (F){
  p1 = link(results$tau[,1])
  p201 = link(results$tau[,201])
  p366 = link(results$tau[,366])
  plot(p1, type="l")
  plot(p201, type="l")
  plot(p366, type="l")
  
  acf(results$tau[,1])
  acf(results$tau[,201])
  acf(results$tau[,365])
  
  hist(p1, nclass=100, prob=T)
  hist(p201, nclass=100, prob=T)
  hist(p366, nclass=100, prob=T)
  
  plot(results$sigma, type="l")
}
par(mfrow=c(1,1))
# gibbs(tau, Q)
