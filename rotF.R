dGetUij <- function(i, j, Kd, A, U, lamdaij) {
  return(sum(Kd[, j] * A[, i]))
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
  newdU <- newdU - new.lambda * ((2*diag(1,2)*new.U)+diag(1,2)[2:1, ])
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
new.U <- eigen(U)$vectors #U
K <- A %*% new.U
Kd <- 4 * K * K * K
new.lambda <- matrix(10:13, 2, 2)
mu.lambda <- 0.01
mu.U <- 0.01
for (iter in 1:2000) {
  new.U <- getU(new.U, mu.U, new.lambda)
  new.U <- eigen(U)$vectors
  new.lambda <- getLambda(new.U, new.lambda, mu.lambda)
  temp <- A %*% new.U
  print(temp)
}

A.rot <- A %*% new.U
plot.new()
lines(x=c(0, A.rot[5, 1]), y=c(0, A.rot[5, 2]), col = "red", xlim)
A.rot

