library("sn"); data(ais)
ais <- ais[,-c(1,2)]
head(ais)
X <- data.matrix(ais[,-8])
head(X)
lm(Bfat~.-1, data=ais)
solve(t(X)%*%X)%*%t(X)%*%ais$Bfat

# Ridge Regression:
alpha=.001
solve(t(X)%*%X+diag(alpha,ncol(X)))%*%t(X)%*%ais$Bfat

# Lasso
library(glmnet)
lam <- .001
est <- glmnet(X, ais$Bfat, lambda = lam)
est$beta
