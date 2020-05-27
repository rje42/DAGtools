A <- matrix(c(0,1,1,0,
              0,0,1,1,
              0,0,0,1,
              0,0,0,0), 4, 4,  byrow = TRUE)

set.seed(129)
Sigma <- genSamp(A, rep(0,4), rWishart(1,10,Sigma=diag(4))[,,1], n=10)

n <- 1e6
dat <- matrix(rnorm(4*n), ncol=4) %*% chol(Sigma)
# S <- cov(dat)

tmp <- causalEffects(A, Sigma)
# print(Sigma)
# tmp$direct_normal

## Look at effect of 1 on 2
rho1 <- cor(dat[,1],dat[,2])
stat1 <- atanh(rho1) - atanh(tmp$direct_normal[1,2])

## Look at effect of 2 on 4
S <- intervene(cov(dat), 2, which(A[,2] > 0))   # not a controlled effect
# rho2 <- (S / sqrt(outer(diag(S),diag(S))))[2,4]    # !!!
rho2 <- cor(lm(dat[,4]~dat[,c(3)])$resid, lm(dat[,2] ~ dat[,1])$resid)
stat2 <- atanh(rho2) - atanh(tmp$direct_normal[2,4])

# print(c(rho1, stat1,  rho2, stat2))

test_that("stats are small", {
  expect_true(abs(stat1) < 4/sqrt(n-3))
  expect_true(abs(stat2) < 4/sqrt(n-4))
})


#### different topology

A <- matrix(c(0,0,1,0,
              1,0,1,1,
              0,0,0,1,
              0,0,0,0), 4, 4,  byrow = TRUE)

set.seed(130)
Sigma <- genSamp(A, rep(0,4), rWishart(1,10,Sigma=diag(4))[,,1], n=100)

n <- 1e6
dat <- matrix(rnorm(4*n), ncol=4) %*% chol(Sigma)
# S <- cov(dat)

tmp <- causalEffects(A, Sigma)
# print(Sigma)
# tmp$direct_normal

## Look at effect of 1 on 2
rho1 <- cor(dat[,1],dat[,2])
stat1 <- atanh(rho1) - atanh(tmp$direct_normal[2,1])

## Look at effect of 2 on 4
S <- intervene(cov(dat), 3, which(A[,3] > 0))   # not a controlled effect
rho2 <- # (S / sqrt(outer(diag(S),diag(S))))[3,4]    # !!!
cor(lm(dat[,4]~dat[,2])$resid, lm(dat[,3] ~ dat[,1]+dat[,2])$resid)
stat2 <- atanh(rho2) - atanh(tmp$direct_normal[3,4])

# print(c(rho1, stat1,  rho2, stat2))

test_that("stats are small", {
  expect_true(abs(stat1) < 4/sqrt(n-3))
  expect_true(abs(stat2) < 10/sqrt(n-4))
})

