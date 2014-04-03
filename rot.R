dGetUij <- function(i, j, Kd, A, U, lamdaij) {
  return(sum(Kd[, j] * A[, i]) - lamdaij[i, j]*U[j, j])
}
dGetlamdaij <- function(i, j, U) {
  U.tilda <- U %*% t(U)
  return ( (i == j) - U.tilda[i, j])
}
getU <- function(new.U, mu.U, new.lambda) {
  newdU <- 0*new.U
  nc <- dim(U)[1]
  nr <- dim(U)[2]
  for (j in 1:nc) {
    for (i in 1:nr){
      newdU[i, j] <- dGetUij(i, j, Kd, A, new.U, new.lambda)
    }
  }
  new.U <- new.U + mu.U * newdU
  new.U
}

getLambda <- function(new.U, new.lambda, mu.lambda) {
  dlambda <- 0*new.lambda
  nc <- dim(new.lambda)[1]
  nr <- dim(new.lambda)[2]
  for (j in 1:nc) {
    for (i in 1:nr){
      dlambda[i, j] <- dGetlamdaij(i, j, new.U)
    }
  }
  new.lambda <- new.lambda + mu.lambda * dlambda
  new.lambda
}

############################################
arr.new <- c(-0.9511, -1.6435, 2.3655, -2.9154, -3.7010, 0.9511, -1.6435, 2.3655, -2.9154, 3.7010)
A <- matrix(arr.new, ncol = 2, nrow = 5)
U <- matrix(1:4,2,2)
U <- U * 0.1
K <- A %*% U
Kd <- 4 * K * K * K
new.U <- eigen(U)$vectors #U
new.lambda <- 0.01 * matrix(10:13, 2, 2)
mu.lambda <- 0.001
mu.U <- 0.001
for (iter in 1:200) {
  new.U <- getU(new.U, mu.U, new.lambda)
  new.U <- eigen(U)$vectors
  new.lambda <- getLambda(new.U, new.lambda, mu.lambda)
  temp <- A %*% new.U
  print(temp)
}


