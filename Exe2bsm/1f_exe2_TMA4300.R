# Libs ----
library(MASS)
library(Matrix)

# 1e) fncs ----
link = function(tau){
  return(exp(tau)/(1+exp(tau)))
}

logbin = function(n, y, tau){
  # Remake and use this
  return(y*log(1+exp(-tau)) - (n-y)*log(1+exp(tau)))
}

acceptRatio = function(n, y, tauProp, tau){
  # Confirmed faster. N=1000: 4.7 vs 3.81
  return(exp(y*(tauProp - tau) + n*log((1+exp(tau))/(1+exp(tauProp)))))
}


mhRW = function(tau, sigma, yt, t, normVec, QaaInv, Qab){
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
  # prop_tau = normVec[t]*sigma + mu_ab
  prop_tau = normVec*sigma + mu_ab
  n = ifelse(t==60,10,39)
  ratio = acceptRatio(n, yt, prop_tau, tau[t])
  if (runif(1) < min(c(1,ratio))){
    return(list(tau=prop_tau, accepted=1))
  }
  else{return(list(tau=tau[t], accepted=0))}
}




mcmcRW = function(N, dt, M = 10, sigma0=0.1){
  # Changed for block
  # Allocate memory
  Ttot = 366
  tau = matrix(NA, nrow=N, ncol = Ttot)
  sigma = numeric(length = N)
  tau_i = numeric(length = Ttot)
  normMat = matrix(rep(rnorm(Ttot), N), nrow = N, ncol=Ttot)
  normVec = normMat[1,]
  
  # Find init vals
  tau[1,] = rnorm(Ttot) # init tau drawn from normal distr.
  # tau[1,] = runif(366, -100, 100) # init tau drawn from uniform distr.
  sigma[1] = sigma0
  
  # Make Q matrices (all inversed Qaa and Qab)
  Qs = QpreComp(Ttot, M)
  Qab = Qs$Qab
  n.sets = Qs$n.sets
  
  
  
  
  # Run mcmc for N iterations
  accepted = 0
  for (i in 2:N){
    tau_i = tau[i-1,]
    sigma_i = sigma[i-1]
    "Add first step using Qaa1"
    rtemp= mhRW(tau_i, sqrt(sigma_i), dt$n.rain[t], t, normVec[t])
    for (t in 1:n.sets){
      "Loop through mid steps using Qaa2"
      rtemp= mhRW(tau_i, sqrt(sigma_i), dt$n.rain[t], t, normVec[t])
      tau[i,t] = rtemp$tau
      accepted = accepted + rtemp$accepted
    }
    "Add last step using Qaa3"
    normVec = normMat[i,]
    # Squared diff. of tau vec.
    tQt = sum((tau[i,-366] - tau[i,-1])^2) # this sim tau vals.
    # tQt = sum((tau[i-1,-366] - tau[i-2,-1])^2) # prev sim tau vals.
    
    # Gibbs step (Draw from IG)
    sigma[i] = 1/rgamma(1, 2 + (366-1)/2, 0.05 + 0.5*tQt) # Gibbs inline
    
    if (i%%(N/10)==0){
      print(i/N*100)
      print(accepted/(i*366))
    }
  }
  return(list(tau=tau, sigma=sigma))
}

# 1f) fncs ----
QpreComp = function(Ttot, M){
  # Precomputes all Q matrices needed for simulation
  # Make sparse Q and Qaa1, Qaa2
  Q <- bandSparse(t, t, #dimensions
                  (-1):1, #band, diagonal is number 0
                  list(rep(-1, Ttot-1), # Tridiag values
                       rep(2, Ttot), 
                       rep(-1, Ttot-1)))
  Q[1,1] = 1 
  Q[t,t] = 1
  
  Qaa1 = Q[1:M,1:M]
  Qaa2 = Q[2:(M+1), 2:(M+1)]
  
  QabList = list(Qab1 = Q[1:M, -(1:M)]) # Qab to be used with Q1inv
  
  # Make intervals and Q3
  n.sets = floor(Ttot/M)
  mod = Ttot %% M
  intervals = matrix(1:(Ttot - mod), ncol = M, byrow=T)
  if (mod){
    nTot = n.sets + 1
    intervalsRes = (Ttot-mod+1):Ttot
    Qaa3 = Q[intervalsRes, intervalsRes]
    for (i in 2:n.sets){
      QabList[i] = Q[intervals[i,], -intervals[i,]]
    }
    QabList[n.sets+1] = Q[intervalsRes, -intervalsRes]
  } else {
    nTot = n.sets
    intervalsRes = NA
    Qaa3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
    for (i in 2:n.sets){
      QabList[i] = Q[intervals[i,], -intervals[i,]]
    }
  }
  
  names(QabList) = sprintf("Qab%d", 1:(n.sets + 1))
  
  return(list(Qaa1inv = solve(Q1), Qaa2inv = solve(Q2), Qaa3inv = solve(Q3),
              Qab = QabList,
              n.sets = nTot))
}
