library(Matrix)

# setwd('./Exe2bsm')
load('./rain.rda')

k = 10
l = list(c(10,2,4),c(23,1))
l
names(l) = c(paste0("tau", k), sigma)

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

Ttot = 366
M = 10

Qs = Qprecomp(Ttot, M)
Qs2 = Qprecomp2(Ttot, M)
all.equal((Qs$Qmult$Qmult1),-Qs2$S[[1]])
all.equal((Qs$Qmult$Qmult2),-Qs2$S[[2]])
all.equal((Qs$Qmult$Qmult3),-Qs2$S[[3]])
all.equal(Qs$Qmult$Qmult37, -Qs2$S[[37]])

all.equal(Qs$Qaa1inv, Qs2$QaaInv[[1]])
all.equal(Qs$Qaa2inv, Qs2$QaaInv[[2]])
all.equal(Qs$Qaa2inv, Qs2$QaaInv[[3]])
all.equal(Qs$Qaa3inv, Qs2$QaaInv[[37]])


# Prev fncs ----
# Qprecomp = function(Ttot, M){
#   # Precomputes all Q matrices needed for simulation
#   ## Make I and Q3
#   n.sets = floor(Ttot/M)
#   mod = Ttot %% M
#   I = matrix(1:(Ttot - mod), ncol = M, byrow=T)
#   ## Make sparse Q and Qaa and Qab
#   # Q <- bandSparse(Ttot, Ttot, #dimensions
#   #                 (-1):1, #band, diagonal is number 0
#   #                 list(rep(-1, Ttot-1), # Tridiag values
#   #                      rep(2, Ttot),
#   #                      rep(-1, Ttot-1)))
#   
#   # Make non sparse Q
#   Q = triDiag(diagonal = 2, upper = -1, lower = -1, nrow = Ttot)
#   
#   Q[1,1] = 1 
#   Q[Ttot,Ttot] = 1
#   
#   Qaa1 = Q[1:M,1:M]
#   Qaa2 = Q[2:(M+1), 2:(M+1)]
#   
#   Qaa1inv = solve(Qaa1)
#   Qaa2inv = solve(Qaa2)
#   
#   # Qab = list(Qab1 = Q[1:M, -(1:M)]) # Qab to be used with Q1inv
#   
#   
#   Qmult = list(QaaInvQab1 = -Qaa1inv %*% Qab$Qab1)
#   
#   
#   if (mod){
#     nTot = n.sets + 1
#     Qab = matrix(NA, nrow = nTot, ncol=Ttot-M)
#     Ires = (Ttot-mod+1):Ttot
#     Qaa3 = Q[Ires, Ires]
#     Qaa3inv = solve(Qaa3)
#     for (i in 2:n.sets){
#       # Qab[i] = Q[I[i,], -I[i,]]
#       Qab = append(Qab,Q[I[i,], -I[i,]])
#       # Qmult[i] = -Qaa2inv %*% Qab[[i]]
#       Qmult = append(Qmult, -Qaa2inv %*% Qab[[i]]) # sparse
#       # Qmult = append(Qmult, -Qaa2inv %*% Qab[i]) # not sparse
#     }
#     # Qab[n.sets+1] = Q[Ires, -Ires]
#     Qab = append(Qab, Q[Ires, -Ires])
#     # Qmult[n.sets + 1] = -Qaa3inv %*% Qab[[n.sets + 1]]
#     Qmult = append(Qmult, -Qaa3inv %*% Qab[[n.sets + 1]]) # sparse
#     # Qmult = append(Qmult, -Qaa3inv %*% Qab[n.sets + 1])
#   } else {
#     nTot = n.sets
#     Qab = matrix(NA, nrow = nTot, ncol=Ttot-M)
#     Qaa3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
#     Qaa3inv = solve(Qaa3)
#     for (i in 2:(n.sets-1)){
#       # Qab[i] = Q[I[i,], -I[i,]]
#       Qab = append(Qab,Q[I[i,], -I[i,]])
#       # Qmult[i] = -Qaa2inv %*% Qab[[i]]
#       Qmult = append(Qmult, -Qaa2inv %*% Qab[[i]]) # sparse
#       # Qmult = append(Qmult, -Qaa2inv %*% Qab[i]) # not sparse
#     }
#     Qab = append(Qab,Q[I[nTot,], -I[nTot,]])
#     Ires = Qab[nTot]
#     Qmult = append(Qmult, -Qaa3inv %*% Qab[[nTot]]) # sparse
#     # Qmult = append(Qmult, -Qaa3inv %*% Qab[nTot]) # not sparse
#   }
#   
#   names(Qab) = sprintf("Qab%d", 1:(n.sets + 1))
#   names(Qmult) = sprintf("Qmult%d", 1:(n.sets + 1))
#   
#   return(list(Qaa1inv = solve(Qaa1), Qaa2inv = solve(Qaa2), Qaa3inv = solve(Qaa3),
#               Qab = Qab, Qmult = Qmult,
#               n.sets = nTot, I = I, Ires = Ires))
# }