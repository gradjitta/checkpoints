
genDataPCA <- function(x, y, rot) {
  mat <-matrix(1000, 1000, 2)
  mat[, 1] <- x
  mat[, 2] <- y
  data <- mat %*% rot
  plot(data[, 1], data[, 2], col = rgb(0, 0, 0, 0.4))
  data
}
pca <- function(data, k) {
  #covariance.mat <- t(data) %*% data
  covariance.mat <- cov(data)
  eig.vecs <- eigen(covariance.mat)$vectors
  eig.vals <- eigen(covariance.mat)$values
  k.components <- eig.vecs[, 1:k]
  k.components
}
covariance.mat <- cov(data)
x <- rnorm(1000, mean = 0, sd = 0.5)
y <- rnorm(1000, mean = 0, sd = 0.5)
rot <-matrix(c(1,2,2,1), 2, 2)
data <- genDataPCA(x, y, rot)
pc_2 <- pca(data, 2)
eig.vecs <- eigen(covariance.mat)$vectors
eig.vals <- eigen(covariance.mat)$values
evecs1 <- pc_2[, 1] * eig.vals[1]
evecs2 <- pc_2[, 2] * eig.vals[2]
lines(x=c(0, evecs1[1]), y=c(0, evecs1[2]), col = "red")
lines(x=c(0, -evecs1[1]), y=c(0, -evecs1[2]), col = "red")
lines(x=c(0, evecs2[1]), y=c(0, evecs2[2]), col = "red")
lines(x=c(0, -evecs2[1]), y=c(0, -evecs2[2]), col = "red")



