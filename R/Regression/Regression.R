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

##Example 2
library(ISLR2)
data("NCI60"); ind <- 1:64
y <- NCI60$data[,6830]; X <- NCI60$data[,ind]
est <- solve(t(X)%*%X)%*%t(X)%*%y
cbind(est, lm(y~X-1)$coeff)

# Ridge Regression:
alpha=0.1
r.est <- solve(t(X)%*%X+diag(alpha,ncol(X)))%*%t(X)%*%y
cbind(est,r.est)

# Lasso
lam <- .1
l.est <- glmnet(X, y, lambda = lam)
l.est$beta
