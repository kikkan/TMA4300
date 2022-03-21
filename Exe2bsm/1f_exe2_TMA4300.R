# Libs ----
library(MASS)
library(Matrix)
library(SoDA)

# 1e) fncs ----
link = function(tau){
  return(exp(tau)/(1+exp(tau)))
}


acceptRatio = function(n, y, tauProp, tau){
  return(exp(y*(tauProp - tau) + n*log((1+exp(tau))/(1+exp(tauProp)))))
}

acceptRatioBlock = function(I, tauProp, tauPrev){
  first = sum(rain$n.rain[I]*(tauProp - tauPrev))
  second = sum(rain$n.years[I]*(log(1+ exp(tauPrev)) - log(1+exp(tauProp))))
  return(exp(first + second))
}


mhRW = function(tau, y, n, tvec, Sigma, QmultI){
  mu_ab = QmultI %*% matrix(tau[-tvec], ncol=1)
  prop_tau = mvrnorm(n=1, mu_ab[,1], Sigma=Sigma)
  
  # ratio = 0
  # for (t in 1:length(n)){
  #   ratio = ratio + acceptRatio(n[t], y[t], prop_tau[t], tau[tvec[t]])
  # }
  # if (runif(1) < min(c(1,ratio))){
  #   return(list(tau=prop_tau, accepted=1))
  # }
  # else{return(list(tau=tau[tvec], accepted=0))}
  
  ratio = acceptRatioBlock(tvec, prop_tau, tau[tvec])
  if (runif(1) < min(c(1,ratio))){
    return(list(tau=prop_tau, accepted=length(tvec)))
  }
  else{return(list(tau=tau[tvec], accepted=0))}
}




mcmcBlock = function(N, dt, M = 10, sigma0=0.1){
  # Changed for block
  # Allocate memory
  Ttot = 366
  tau = matrix(NA, nrow = N, ncol = Ttot) 
  sigma = numeric(length = N)
  tauOld = numeric(length = Ttot) # Remove?
  
  # Find init vals
  tau[1,] = rnorm(Ttot) # init tau drawn from normal distr.
  # tau[1,] = runif(366, -100, 100) # init tau drawn from uniform distr.
  sigma[1] = sigma0
  
  # Make Q matrices (all inversed Qaa and Qab)
  Qs = Qprecomp(Ttot, M)
  Qmult = Qs$Qmult # Precomp. of - inv Qaa times Qab
  n.sets = Qs$n.sets
  I = Qs$I
  Ires = Qs$Ires
  
  # Run mcmc for N iterations
  accepted = 0
  for (i in 2:N){
    # First step using Qaa1inv
    rtemp= mhRW(tau[i-1,], y=dt$n.rain[I[,1]], n=dt$n.years[I[,1]], 
                tvec=I[,1], sigma[i-1]*Qs$Qaa1inv, Qmult[[1]])
    tau[i, I[,1]] = rtemp$tau
    accepted = accepted + rtemp$accepted
    
    Sigma2 = sigma[i-1]*Qs$Qaa2inv # covar mat for 1<t<366
    for (t in 2:(n.sets-1)){
      rtemp= mhRW(tau=tau[i-1,], y=dt$n.rain[I[,t]], n=dt$n.years[I[,t]], 
                  tvec=I[,t], Sigma=Sigma2, QmultI=Qmult[[t]])
      tau[i, I[,t]] = rtemp$tau
      accepted = accepted + rtemp$accepted
    }
    
    rtemp= mhRW(tau[i-1,], dt$n.rain[Ires], dt$n.years[Ires], 
                Ires, sigma[i-1]*Qs$Qaa3inv, Qmult[[n.sets]])
    tau[i, Ires] = rtemp$tau
    accepted = accepted + rtemp$accepted
    
    # tQt = sum((tau[i,-366] - tau[i,-1])^2) 
    tQt = sum((diff(tau[i,]))^2) 
    sigma[i] = 1/rgamma(1, shape=2 + (Ttot-1)/2, rate=0.05 + 0.5*tQt) # Gibbs inline
    
    if (i%%(N/10)==0){
      print(i/N*100)
      print(accepted/(i*366))
    }
  }
  return(list(tau=tau, sigma=sigma))
}



# 1f) fncs ----
Qprecomp = function(Ttot, M){
  # Precomputes all Q matrices needed for simulation
  ## Make I and Q3
  n.sets = ceiling(Ttot/M)
  mod = Ttot %% M
  I = matrix(1:(M*(n.sets-1)), ncol = n.sets-1, byrow=F)
  Ires = (I[M,n.sets-1]+1):Ttot
  ## Make sparse Q and Qaa and Qab
  # Q <- bandSparse(Ttot, Ttot, #dimensions
  #                 (-1):1, #band, diagonal is number 0
  #                 list(rep(-1, Ttot-1), # Tridiag values
  #                      rep(2, Ttot),
  #                      rep(-1, Ttot-1)))
  
  # Make non sparse Q
  Q = triDiag(diagonal = 2, upper = -1, lower = -1, nrow = Ttot)
  
  Q[1,1] = 1 
  Q[Ttot,Ttot] = 1
  
  Qaa1 = Q[1:M,1:M]
  Qaa2 = Q[2:(M+1), 2:(M+1)]
  Qaa3 = Q[Ires, Ires]
  
  Qaa1inv = solve(Qaa1)
  Qaa2inv = solve(Qaa2)
  Qaa3inv = solve(Qaa3)
  
  Qmult = list(Qmult1 = -Qaa1inv %*% Q[1:M, -(1:M)])
  for (i in 2:(n.sets-1)){
    Qmult = append(Qmult, list(-Qaa2inv %*% Q[I[,i], -I[,i]]))
  }
  Qmult = append(Qmult, list(-Qaa3inv %*% Q[Ires, -Ires]))
  names(Qmult) = sprintf("Qmult%d", 1:(n.sets))
  
  return(list(Qaa1inv = Qaa1inv, Qaa2inv = Qaa2inv, Qaa3inv = Qaa3inv,
              Qmult = Qmult,
              n.sets = n.sets, I = I, Ires = Ires))
}


if (F){
  load('./rain.rda')
  setwd('./Exe2bsm/')
  options(error = recover)
  options(error = NULL)
}

M = 10
N = 10000

ptm = proc.time()
# Run ----
set.seed(321)
results = mcmcBlock(N, rain, M)
proc.time() - ptm

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
  
  par(mfrow=c(1,1))
  plot(results$sigma, type="l")
}

par(mfrow=c(3,1))
plot(link(results$tau[1,]), type = "l")
plot(link(results$tau[N/2,]), type = "l")
plot(link(results$tau[N,]), type = "l")
