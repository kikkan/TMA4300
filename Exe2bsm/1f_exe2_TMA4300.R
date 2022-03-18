# Libs ----
library(MASS)
library(Matrix)

# 1e) fncs ----
link = function(tau){
  return(exp(tau)/(1+exp(tau)))
}


acceptRatio = function(n, y, tauProp, tau){
  return(exp(y*(tauProp - tau) + n*log((1+exp(tau))/(1+exp(tauProp)))))
}


mhRW = function(tau, y, n, tvec, Sigma, QmultI){
  mu_ab = QmultI %*% matrix(tau[-tvec], ncol=1)
  prop_tau = mvrnorm(n=1, mu_ab[,1], Sigma=Sigma)
  ratio = 0
  for (t in 1:length(n)){
    ratio = ratio + acceptRatio(n[t], y[t], prop_tau[t], tau[tvec[t]])
  }
  if (runif(1) < min(c(1,ratio))){
    return(list(tau=prop_tau, accepted=1))
  }
  else{return(list(tau=tau[tvec], accepted=0))}
}




mcmcBlock = function(N, dt, M = 10, sigma0=0.1){
  # Changed for block
  # Allocate memory
  Ttot = 366
  tau = matrix(NA, nrow=N, ncol = Ttot) 
  sigma = numeric(length = N)
  tauOld = numeric(length = Ttot)
  normMat = matrix(rep(rnorm(Ttot), N), nrow = N, ncol=Ttot)
  normVec = normMat[1,]
  
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
    tauOld = tau[i-1,]
    sigma_i = sigma[i-1]
    # First step using Qaa1inv
    # browser()
    Sigma2 = sigma_i*Qs$Qaa2inv # covar mat for 1<t<366
    rtemp= mhRW(tauOld, y=dt$n.rain[I[1,]], n=dt$n.years[I[1,]], 
                tvec=I[1,], sigma_i*Qs$Qaa1inv, Qmult[[1]])
    tau[i, I[1,]] = rtemp$tau
    accepted = accepted + rtemp$accepted
    
    for (t in 2:(n.sets-1)){
      rtemp= mhRW(tau=tauOld, y=dt$n.rain[I[t,]], n=dt$n.years[I[t,]], 
                  tvec=I[t,], Sigma=Sigma2, QmultI=Qmult[[t]])
      tau[i, I[t,]] = rtemp$tau
      accepted = accepted + rtemp$accepted
    }
    
    rtemp= mhRW(tauOld, dt$n.rain[Ires], dt$n.years[Ires], 
                Ires, sigma_i*Qs$Qaa3inv, Qmult[[n.sets]])
    tau[i, Ires] = rtemp$tau
    accepted = accepted + rtemp$accepted
    
    # tQt = sum((tau[i,-366] - tau[i,-1])^2) 
    tQt = sum((diff(tau[i,]))^2) 
    sigma[i] = 1/rgamma(1, shape=2 + (366-1)/2, rate=0.05 + 0.5*tQt) # Gibbs inline
    
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
  n.sets = floor(Ttot/M)
  mod = Ttot %% M
  I = matrix(1:(Ttot - mod), ncol = M, byrow=T)
  ## Make sparse Q and Qaa and Qab
  Q <- bandSparse(Ttot, Ttot, #dimensions
                  (-1):1, #band, diagonal is number 0
                  list(rep(-1, Ttot-1), # Tridiag values
                       rep(2, Ttot),
                       rep(-1, Ttot-1)))
  
  # Make non sparse Q
  # Q=matrix(0, nrow = 366, ncol = 366)
  # diag(Q)=2
  # Q[c(1, length(Q))]=1
  # Q[abs(row(Q) - col(Q)) == 1] <- -1
  
  Q[1,1] = 1 
  Q[Ttot,Ttot] = 1
  
  Qaa1 = Q[1:M,1:M]
  Qaa2 = Q[2:(M+1), 2:(M+1)]
  
  Qaa1inv = solve(Qaa1)
  Qaa2inv = solve(Qaa2)
  
  Qab = list(Qab1 = Q[1:M, -(1:M)]) # Qab to be used with Q1inv
  
  Qmult = list(QaaInvQab1 = -Qaa1inv %*% Qab$Qab1)
  
  
  if (mod){
    nTot = n.sets + 1
    Ires = (Ttot-mod+1):Ttot
    Qaa3 = Q[Ires, Ires]
    Qaa3inv = solve(Qaa3)
    for (i in 2:n.sets){
      # Qab[i] = Q[I[i,], -I[i,]]
      Qab = append(Qab,Q[I[i,], -I[i,]])
      # Qmult[i] = -Qaa2inv %*% Qab[[i]]
      Qmult = append(Qmult, -Qaa2inv %*% Qab[[i]]) # sparse
      # Qmult = append(Qmult, -Qaa2inv %*% Qab[i]) # not sparse
    }
    # Qab[n.sets+1] = Q[Ires, -Ires]
    Qab = append(Qab, Q[Ires, -Ires])
    # Qmult[n.sets + 1] = -Qaa3inv %*% Qab[[n.sets + 1]]
    Qmult = append(Qmult, -Qaa3inv %*% Qab[[n.sets + 1]]) # sparse
    # Qmult = append(Qmult, -Qaa3inv %*% Qab[n.sets + 1])
  } else {
    nTot = n.sets
    Qaa3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
    Qaa3inv = solve(Qaa3)
    for (i in 2:(n.sets-1)){
      # Qab[i] = Q[I[i,], -I[i,]]
      Qab = append(Qab,Q[I[i,], -I[i,]])
      # Qmult[i] = -Qaa2inv %*% Qab[[i]]
      Qmult = append(Qmult, -Qaa2inv %*% Qab[[i]]) # sparse
      # Qmult = append(Qmult, -Qaa2inv %*% Qab[i]) # not sparse
    }
    Qab = append(Qab,Q[I[nTot,], -I[nTot,]])
    Ires = Qab[nTot]
    Qmult = append(Qmult, -Qaa3inv %*% Qab[[nTot]]) # sparse
    # Qmult = append(Qmult, -Qaa3inv %*% Qab[nTot]) # not sparse
  }
  
  names(Qab) = sprintf("Qab%d", 1:(n.sets + 1))
  names(Qmult) = sprintf("Qmult%d", 1:(n.sets + 1))
  
  return(list(Qaa1inv = solve(Qaa1), Qaa2inv = solve(Qaa2), Qaa3inv = solve(Qaa3),
              Qab = Qab, Qmult = Qmult,
              n.sets = nTot, I = I, Ires = Ires))
}


if (F){
  load('./rain.rda')
  setwd('./Exe2bsm/')
  options(error = recover)
  options(error = NULL)
}

M = 10
N = 1000

set.seed(321)
# Run ----
ptm = proc.time()
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
