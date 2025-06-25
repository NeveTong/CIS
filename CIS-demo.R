library(Rcpp)
sourceCpp("CI.cpp")
library(Compositional)
n <- 100 
p <- 4   
s <- 1
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- matrix(rnorm(n), nrow = n, ncol = 1)
Z <- matrix(rnorm(n * s), nrow = n, ncol = s)
h <- density(Z)$bw         # kernel bandwidth, for one-dimensional Z
# h <- mkde.tune(Z)$hopt   # for two-dimensional Z

CIS(X, Y, Z, h)          # CIS marginal utility measure
