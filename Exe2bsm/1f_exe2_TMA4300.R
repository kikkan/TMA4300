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


mhRW = function(tau, sigma, yt, t, normVec=NA){
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
  
  # Make Q matrices
  Q=matrix(0, nrow = Ttot, ncol = Ttot)
  diag(Q)=2
  Q[c(1, length(Q))]=1
  Q[abs(row(Q) - col(Q)) == 1] <- -1
  
  # Make sparse Q
  Q <- bandSparse(Ttot, Ttot, #dimensions
                  (-1):1, #band, diagonal is number 0
                  list(rep(-1, Ttot-1), # Tridiag values
                       rep(2, Ttot), 
                       rep(-1, Ttot-1)))
  Q[1,1] = 1 
  Q[t,t] = 1
  
  Qaa1 = Q[1:M,1:M]
  Qaa2 = Q[2:(M+1), 2:(M+1)]
  Qaa3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
  
  # setup all indexes to not have if statement 50000
  
  
  # Run mcmc for N iterations
  accepted = 0
  for (i in 2:N){
    tau_i = tau[i-1,]
    sigma_i = sigma[i-1]
    for (t in 1:366){
      # rtemp= mhRW(tau_i, sigma_i, dt$n.rain[t], t)
      rtemp= mhRW(tau_i, sqrt(sigma_i), dt$n.rain[t], t, normVec[t])
      tau[i,t] = rtemp$tau
      accepted = accepted + rtemp$accepted
    }
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
