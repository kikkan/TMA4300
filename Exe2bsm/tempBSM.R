library(Matrix)


t = 30
M = 5

Q=matrix(0, nrow = t, ncol = t)
diag(Q)=2
Q[c(1, length(Q))]=1
Q[abs(row(Q) - col(Q)) == 1] <- -1

Q1 = Q[1:M,1:M]
Q2 = Q[2:(M+1), 2:(M+1)]
Q3 = Q[(t-M):t, (t-M):t]

Cholesky(Q1)


Q <- bandSparse(t, t, #dimensions
                (-1):1, #band, diagonal is number 0
                list(rep(-1, t-1), # Tridiag values
                     rep(2, t), 
                     rep(-1, t-1)))
Q[1,1] = 1 
Q[t,t] = 1

Q1 = Q[1:M,1:M]
Q2 = Q[2:(M+1), 2:(M+1)]
Q3 = Q[(t-M):t, (t-M):t]

Q1
Cholesky(Q1)
solve(Q2)
solve(Q3)

c = seq(1:20)
ab = c(1,2,3,4,5)
for (i in 1:4){
  print(c[ab])
  ab = ab + 5
}