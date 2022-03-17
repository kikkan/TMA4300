library(Matrix)


Ttot = 30
M = 7

# Make sparse Q and Qaa1, Qaa2
Q <- bandSparse(t, t, #dimensions
                (-1):1, #band, diagonal is number 0
                list(rep(-1, Ttot-1), # Tridiag values
                     rep(2, Ttot), 
                     rep(-1, Ttot-1)))
Q[1,1] = 1 
Q[t,t] = 1

Q1 = Q[1:M,1:M]
Q2 = Q[2:(M+1), 2:(M+1)]

QabList = list(Qab1 = Q[1:M, -(1:M)])



# Make intervals and Q3
n.sets = floor(Ttot/M)
mod = Ttot %% M
intervals = matrix(1:(Ttot - mod), ncol = M, byrow=T)
if (mod){
  intervalsRes = (Ttot-mod+1):Ttot
  Q3 = Q[intervalsRes, intervalsRes]
  for (i in 2:n.sets){
    QabList[i] = Q[intervals[i,], -intervals[i,]]
  }
  QabList[n.sets+1] = Q[intervalsRes, -intervalsRes]
} else {
  intervalsRes = NA
  Q3 = Q[(Ttot-M):Ttot, (Ttot-M):Ttot]
  for (i in 2:n.sets){
    QabList[i] = Q[intervals[i,], -intervals[i,]]
  }
}

names(QabList) = sprintf("Qab[%d]", 1:(n.sets + 1))




for (q in QabList){
  print(q)
}

