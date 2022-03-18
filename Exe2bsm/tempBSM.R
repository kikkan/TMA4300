library(Matrix)

# setwd('./Exe2bsm')
load('./rain.rda')


QpreComp = function(Ttot, M){
  # Precomputes all Q matrices needed for simulation
  # Make sparse Q and Qaa1, Qaa2
  Q <- bandSparse(Ttot, Ttot, #dimensions
                  (-1):1, #band, diagonal is number 0
                  list(rep(-1, Ttot-1), # Tridiag values
                       rep(2, Ttot), 
                       rep(-1, Ttot-1)))
  Q[1,1] = 1 
  Q[Ttot,Ttot] = 1
  
  Qaa1 = Q[1:M,1:M]
  Qaa2 = Q[2:(M+1), 2:(M+1)]
  
  Qaa1inv = solve(Qaa1)
  Qaa2inv = solve(Qaa2)
  
  QabList = list(Qab1 = Q[1:M, -(1:M)]) # Qab to be used with Q1inv
  
  QaaInvQabList = list(QaaInvQab1 = -Qaa1inv %*% QabList$Qab1)
  
  # Make intervals and Q3
  n.sets = floor(Ttot/M)
  mod = Ttot %% M
  intervals = matrix(1:(Ttot - mod), ncol = M, byrow=T)
  if (mod){
    nTot = n.sets + 1
    intervalsRes = (Ttot-mod+1):Ttot
    Qaa3 = Q[intervalsRes, intervalsRes]
    Qaa3inv = solve(Qaa3)
    for (i in 2:n.sets){
      QabList[i] = Q[intervals[i,], -intervals[i,]]
      QaaInvQabList[i] = -Qaa2inv %*% QabList[[i]]
    }
    QabList[n.sets+1] = Q[intervalsRes, -intervalsRes]
    QaaInvQabList[n.sets + 1] = -Qaa3inv %*% QabList[[n.sets + 1]]
  } else {
    nTot = n.sets
    Qaa3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
    for (i in 2:n.sets){
      QabList[i] = Q[intervals[i,], -intervals[i,]]
      QaaInvQabList[i] = -Qaa2inv %*% QabList[[i]]
    }
    intervalsRes = QabList[n.sets]
  }
  
  names(QabList) = sprintf("Qab%d", 1:(n.sets + 1))
  names(QaaInvQabList) = sprintf("QaaInvQab%d", 1:(n.sets + 1))
  
  return(list(Qaa1 = Qaa1, Qaa2 = Qaa2, Qaa3 = Qaa3, # Maybe not needed
              Qaa1inv = solve(Qaa1), Qaa2inv = solve(Qaa2), Qaa3inv = solve(Qaa3),
              Qab = QabList, QaaInvQab = QaaInvQabList,
              n.sets = nTot, intervals = intervals, intRes = intervalsRes))
}

Ttot = 30
M = 7

Qs = QpreComp(Ttot, M)
Qab = Qs$Qab
QaaQab = Qs$QaaInvQab



Qs$Qaa1inv %*% Qab$Qab1
QaaQab$QaaInvQab1
Qs$Qaa2inv %*% Qab$Qab2
-QaaQab$QaaInvQab2
Qs$Qaa2inv %*% Qab$Qab3
Qs$Qaa2inv %*% Qab$Qab4
Qs$Qaa3inv %*% Qab$Qab5
Qab$Qab1
Qs$n.sets

# 1f) Acceptance ratio ----
link = function(tau){
  return(exp(tau)/(1+exp(tau)))
}

tau = rnorm(10)
y = rain$n.rain[58:68]
n = rain$n.years[58:68]

dbinom(y, n, link(tau))


# Testing ----
M = 10
Ttot = 366
l = list(1, ceiling(Ttot/10)-2, 1)
Ttot-(ceiling(Ttot/M)-1)*M
Ttot%%M
