
library(nimble)
library(nimbleEcology)


##################################################
############ Capture-Recapture Models ############
##################################################

## number of individuals
N <- 800
## number of observation periods
T <- 8
## time period of first capture
first <- rep(1, N)
## length of observation history from first capture
len <- T - first + 1
## survival probability
phi <- 0.7
## detection probability
pVec <- rep(c(0.2, 0.8), each=N/2)

## simulate z (alive / dead status),
## and y (encounter histories)
set.seed(0)
z <- matrix(NA, nrow=N, ncol=T)
y <- matrix(NA, nrow=N, ncol=T)
for(i in 1:N) {
    z[i, first[i]] <- y[i, first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1, phi*z[i,t-1])
        y[i,t] <- rbinom(1, 1, pVec[i]*z[i,t])
    }
}

## Homogeneous Capture-Recapture Model
code <-  nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        y[i,first[i]:T] ~ dCJS_ss(phi, p, len=len[i])
    }
})
constants <- list(N=N, T=T, first=first, len=len)
data <- list(y=y)
inits <- list(phi=0.5, p=0.5)
Rmodel <- nimbleModel(code, constants, data, inits)

## 2-Group Finite Mixture Capture-Recapture Model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    for(k in 1:K)   p[k] ~ dunif(0, 1)
    one ~ dconstraint(p[1] <= p[2])
    for(i in 1:N) {
        g[i] ~ dcat(pi[1:K])
        y[i,first[i]:T] ~ dCJS_ss(phi, p[g[i]], len=len[i])
    }
})
K <- 2    ## fixed number of groups
constants <- list(N=N, T=T, first=first, len=len, K=K)
data <- list(y=y, one=rep(1,K-1))
inits <- list(phi=0.5, pi=rep(1/K,K), p=rep(0.5,K), g=rep(1,N))
Rmodel <- nimbleModel(code, constants, data, inits)

## 3-Group Finite Mixture Capture-Recapture Model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    for(k in 1:K)   p[k] ~ dunif(0, 1)
    one[1] ~ dconstraint(p[1] <= p[2])
    one[2] ~ dconstraint(p[2] <= p[3])
    for(i in 1:N) {
        g[i] ~ dcat(pi[1:K])
        y[i,first[i]:T] ~ dCJS_ss(phi, p[g[i]], len=len[i])
    }
})
K <- 3    ## fixed number of groups
constants <- list(N=N, T=T, first=first, len=len, K=K)
data <- list(y=y, one=rep(1,K-1))
inits <- list(phi=0.5, pi=rep(1/K,K), p=rep(0.5,K), g=rep(1,N))
Rmodel <- nimbleModel(code, constants, data, inits)

## Non-Parametric Capture-Recapture Model
code <- nimbleCode({
    phi ~ dunif(0, 1)
    alpha ~ dgamma(1, 1)
    xi[1:N] ~ dCRP(conc=alpha, size=N)
    for(i in 1:M)   p[i] ~ dunif(0, 1)
    for(i in 1:N) {
        y[i,first[i]:T] ~ dCJS_ss(phi, p[xi[i]], len=len[i])
    }
})
M <- 100   ## maximum number of subgroups
constants <- list(N=N, T=T, first=first, len=len, M=M)
data <- list(y=y)
inits <- list(phi=0.5, alpha=1, xi=rep(1,N), p=rep(0.5,M))
Rmodel <- nimbleModel(code, constants, data, inits)


##################################################
################ Occupancy Models ################
##################################################

## number of sites
N <- 4000
## number of observation periods
T <- 6
## probability of occupancy
pOcc <- 0.7
## detection probability
pVec <- rep(c(0.2, 0.8), each=N/2)

## simulate z (occupied status),
## and y (encounter histories)
set.seed(0)
z <- rep(NA, N)
y <- matrix(NA, nrow=N, ncol=T)
for(i in 1:N) {
    z[i] <- rbinom(1, size=1, prob=pOcc)
    y[i, 1:T] <- rbinom(T, size=1, prob=z[i]*pVec[i])
}

## Homogeneous Occupancy Model
code <- nimbleCode({
    pOcc ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        y[i,1:T] ~ dOcc_s(pOcc, p, len=T)
    }
})
constants <- list(N=N, T=T)
data <- list(y=y)
inits <- list(pOcc=0.5, p=0.5)
Rmodel <- nimbleModel(code, constants, data, inits)

## 2-Group Finite Mixture Occupancy Model
code <- nimbleCode({
    pOcc ~ dunif(0, 1)
    for(k in 1:K)   p[k] ~ dunif(0, 1)
    one ~ dconstraint(p[1] <= p[2])
    for(i in 1:N) {
        g[i] ~ dcat(pi[1:K])
        y[i,1:T] ~ dOcc_s(pOcc, p[g[i]], len=T)
    }
})
K <- 2    ## fixed number of groups
constants <- list(N=N, T=T, K=K)
data <- list(y=y, one=rep(1,K-1))
inits <- list(pOcc=0.5, pi=rep(1/K,K), p=rep(0.5,K), g=rep(1,N))
Rmodel <- nimbleModel(code, constants, data, inits)

## 3-Group Finite Mixture Occupancy Model
code <- nimbleCode({
    pOcc ~ dunif(0, 1)
    for(k in 1:K)   p[k] ~ dunif(0, 1)
    one[1] ~ dconstraint(p[1] <= p[2])
    one[2] ~ dconstraint(p[2] <= p[3])
    for(i in 1:N) {
        g[i] ~ dcat(pi[1:K])
        y[i,1:T] ~ dOcc_s(pOcc, p[g[i]], len=T)
    }
})
K <- 3    ## fixed number of groups
constants <- list(N=N, T=T, K=K)
data <- list(y=y, one=rep(1,K-1))
inits <- list(pOcc=0.5, pi=rep(1/K,K), p=rep(0.5,K), g=rep(1,N))
Rmodel <- nimbleModel(code, constants, data, inits)

## Non-Parametric Occupancy Model
code <- nimbleCode({
    pOcc ~ dunif(0, 1)
    alpha ~ dgamma(1, 1)
    xi[1:N] ~ dCRP(conc=alpha, size=N)
    for(i in 1:M)   p[i] ~ dunif(0, 1)
    for(i in 1:N) {
        y[i,1:T] ~ dOcc_s(pOcc, p[xi[i]], len=T)
    }
})
M <- 100   ## maximum number of subgroups
constants <- list(N=N, T=T, M=M)
data <- list(y=y)
inits <- list(pOcc=0.5, alpha=1, xi=rep(1,N), p=rep(0.5,M))
Rmodel <- nimbleModel(code, constants, data, inits)


##################################################
############## Fit Model Using MCMC ##############
##################################################

## configure MCMC
conf <- configureMCMC(Rmodel)

## build MCMC
Rmcmc <- buildMCMC(conf)

## compile model and MCMC
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)

set.seed(0)
samplesList <- runMCMC(Cmcmc, niter=10000, nchains=3)
