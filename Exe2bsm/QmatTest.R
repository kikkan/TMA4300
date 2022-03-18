Qprecomp2 <- function(T, M){
  # Blocking is now possible for M>1
  # if(M==1){
  #   return( mcmc.iterative(num.iter, sigma0, tau0) )
  # }
  ## Initialising
  # tau.mat <- matrix(NA, nrow = num.iter, ncol = T)
  # tau.mat[1, ] <- tau0
  # sigma.vec <- rep(NA, num.iter)
  # sigma.vec[1] <- sigma0
  # count <- 0 # Count of accepted tau-samples
  # alpha.vec <- rep(NA, num.iter - 1) # store acc 
  
  
  ## Precomputing, (assuming M<T)
  n.block.total <- ceiling(T/M)
  n.blocks <- list(1, ceiling(T/M)-2, 1)
  M.3 <- T-(ceiling(T/M)-1)*M # mod
  # Q.AA for the three different blocks
  Q.AA <- list( Q[1:M, 1:M], Q[2:(M+1), 2:(M+1)], Q[(T-M.3+1):T, (T-M.3+1):T] )
  # Inverse of Q.AA for all blocks
  Q.AA.inv <- list( solve(Q.AA[[1]]) )
  if(n.blocks[[2]] > 0){
    Q.AA.inv2 <- solve(Q.AA[[2]])
    for(n in 1:n.blocks[[2]]){
      Q.AA.inv <- append( Q.AA.inv, list(Q.AA.inv2) )
    }  
  } 
  Q.AA.inv <- append( Q.AA.inv, list(solve(Q.AA[[3]])) )
  # Q.AB for all blocks
  Q.AB <- list(Q[1:M, (M+1):T])
  if(n.blocks[[2]] > 0){
    for(n in 1:n.blocks[[2]]){
      I <- (M*n + 1):(M*(n+1))
      Q.AB <- append( Q.AB, list(Q[I, (1:T)[-I]]) )
    }   
  }
  Q.AB <- append( Q.AB, list(Q[(T-M.3+1):T, 1:(T-M.3)]) )
  # inv(Q.AA) * Q.AB
  S <- list()
  for(n in 1:n.block.total){
    S <- append( S, list(Q.AA.inv[[n]] %*% Q.AB[[n]]) )
  }
  return(list(S = S, Qab=Q.AB, QaaInv = Q.AA.inv))
}



Q <- bandSparse(Ttot, Ttot, #dimensions
                (-1):1, #band, diagonal is number 0
                list(rep(-1, Ttot-1), # Tridiag values
                     rep(2, Ttot), 
                     rep(-1, Ttot-1)))
Q[1,1] = 1 
Q[Ttot,Ttot] = 1 

Q1 = Qprecomp(366, 10)
Q2 = Qprecomp2(366, 6)

last = length(Q2$QaaInv)

Q1$Qaa1inv == Q2$QaaInv[[1]]
Q1$Qaa2inv == Q2$QaaInv[[2]]
Q1$Qaa2inv == Q2$QaaInv[[3]]
Q1$Qaa3inv == Q2$QaaInv[[last]]
Q2$QaaInv[last]

Q1$Qmult$Qmult2





all.equal(-Q1$Qmult$Qmult1, Q2$S[[1]])
