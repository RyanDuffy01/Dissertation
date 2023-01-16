#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(nimble)
library(igraph)
library(coda)
library(R6)

# This is an example of how to use Nimble and the BUGS language for a simple conjugate analysis 

# Revisit the first example in Section 1.2.2 of your notes. 
# X[i] is the lifetime of laptop i, i=1,...,n. 
# Exponential likelihood given \lambda. X[i]~Exp(\lambda)
# 1/lambda is the average lifetime of a certain laptop brand
# Assume n=20, with \xbar=5
# Gamma prior on \lambda

################################################################
################################################################

# Specify the statistical model
ExpGammaCode <- nimbleCode({
  
 # Specify the likelihood:
  for (i in 1:N){
    x[i] ~ dexp(lambda)
  }
 # to generate values for the average lifetime, given the \lambda samples
  averagelifetime<-1/lambda
  
 # Prior specification:
  lambda ~ dgamma(c,d) # Gamma prior with known parameters c and d
  
})

# Values for some constants in the model
ExpGammaConsts <- list(N = 20, c=0.2, d=0.6) # for a vague prior
# ExpGammaConsts <- list(N = 20, c=10, d=30) # for an informative prior

# The data values
ExpGammaData <- list(x=c(4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6))

# one set of initial values before building the model                 
ExpGammaInits <- list(lambda=1)  # missing data are random variables and need to be initialised too. 
                                 # Nimble can do this by sampling from the priors. 

# to build the model
ExpGamma <- nimbleModel(code = ExpGammaCode, name = "ExpGamma", constants = ExpGammaConsts,
                    data = ExpGammaData, inits= ExpGammaInits)

# To compile the model
CExpGamma <- compileNimble(ExpGamma)

# set up the monitored quantities. Default is all of the random quantities
ExpGammaConf <- configureMCMC(ExpGamma, monitors = c('lambda','averagelifetime'), print = TRUE) 
# build the MCMC algorithm
ExpGammaMCMC <- buildMCMC(ExpGammaConf)
# compile the MCMC chain 
CExpGammaMCMC <- compileNimble(ExpGammaMCMC, project = ExpGamma)

####################################################################################
####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS  #
####################################################################################
set.seed(10)
ExpGammaInits <- list(list(lambda = 1), list(lambda = 10))
posterior <- runMCMC(CExpGammaMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = ExpGammaInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
windows()# or quartz() on a MacOS
plot(combinedchains)
#gelman.diag(combinedchains)
windows();gelman.plot(combinedchains)
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
effectiveSize(combinedchains)
posterior$summary$all.chains
summary(combinedchains)# returns Time-series SE


###################################################################################
##### TO PREDICT THE LIFETIME OF A LAPTOP BOUGHT IN THE FUTURE   ##################
###################################################################################
# Specify the statistical model
ExpGammaCode <- nimbleCode({
  # Specify the likelihood:
  for (i in 1:N){
    x[i] ~ dexp( lambda)
  }
  # to generate values for the average lifetime, given the \lambda samples
  x.new~ dexp( lambda) # Nimble will sample from the posterior predictive distribution
  averagelifetime<-1/lambda
  # Prior specification:
  lambda~dgamma(10,30) # for an informative prior
})
# Values for some constants in the model
ExpGammaConsts <- list(N = 20)
# The data values
ExpGammaData <- list(x=c(4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6))
# one set of initial values before building the model                 
ExpGammaInits <- list(lambda=1, x.new=2)  
# to build the model
ExpGamma <- nimbleModel(code = ExpGammaCode, name = "ExpGamma", constants = ExpGammaConsts,
                        data = ExpGammaData, inits<-ExpGammaInits)
# To compile the model
CExpGamma <- compileNimble(ExpGamma)
# set up the monitored quantities. Default is all of the random quantities
ExpGammaConf <- configureMCMC(ExpGamma, monitors = c('lambda','averagelifetime','x.new'), print = TRUE) 
# build the MCMC algorithm
ExpGammaMCMC <- buildMCMC(ExpGammaConf)
# compile the MCMC chain 
CExpGammaMCMC <- compileNimble(ExpGammaMCMC, project = ExpGamma)
# To obtain the posterior samples
set.seed(10)
ExpGammaInits <- list(list(lambda = 1, x.new=2), list(lambda = 10, x.new=10))
posterior <- runMCMC(CExpGammaMCMC, niter = 100000, thin=1, nburnin=100, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = ExpGammaInits) 
posterior$summary$all.chains

###################################################################################
##### TO PREDICT THE LIFETIME OF A LAPTOP BOUGHT IN THE FUTURE   ##################
##### TREATING THIS LAPTOP AS A MISSING OBSERVATION (IDENTICAL RESULTS AS ABOVE!) #
###################################################################################
# Specify the statistical model
ExpGammaCode <- nimbleCode({
  # Specify the likelihood:
  for (i in 1:N){
    x[i] ~ dexp( lambda)
  }
  # to generate values for the average lifetime, given the \lambda samples
  averagelifetime<-1/lambda
  # Prior specification:
  lambda~dgamma(10,30) # for an informative prior
})
# Values for some constants in the model
ExpGammaConsts <- list(N = 21)
# The data values
ExpGammaData <- list(x=c(4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,NA))
# one set of initial values before building the model                 
ExpGammaInits <- list(lambda=1,x=c(rep(NA,20),3))  
# to build the model
ExpGamma <- nimbleModel(code = ExpGammaCode, name = "ExpGamma", constants = ExpGammaConsts,
                        data = ExpGammaData, inits<-ExpGammaInits)
# To compile the model
CExpGamma <- compileNimble(ExpGamma)
# set up the monitored quantities. Default is all of the random quantities
ExpGammaConf <- configureMCMC(ExpGamma, monitors = c('lambda','averagelifetime','x[21]'), print = TRUE) 
# build the MCMC algorithm
ExpGammaMCMC <- buildMCMC(ExpGammaConf)
# compile the MCMC chain 
CExpGammaMCMC <- compileNimble(ExpGammaMCMC, project = ExpGamma)
# To obtain the posterior samples
set.seed(10)
ExpGammaInits <- list(
      list(lambda = 1, x=c(rep(NA,20),3)),
      list(lambda = 10,x=c(rep(NA,20),9)))
posterior <- runMCMC(CExpGammaMCMC, niter = 1000, thin=1, nburnin=100, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = ExpGammaInits) 
posterior$summary$all.chains
#
########################################################################
####### POSTERIOR SAMPLES AS R MATRICES   ##############################
########################################################################

par(mfrow = c(1, 1))
plot(posterior$samples$chain1[ , "lambda"], type = "l", xlab = "iteration",
     ylab = expression(alpha))
lines(posterior$samples$chain2[ , "lambda"], type = "l", xlab = "iteration", col=2, 
      ylab = expression(alpha))
par(mfrow = c(1, 1))
acf(posterior$samples$chain1[, "lambda"]) # plot autocorrelation of lambda sample - chain 1
par(mfrow = c(1, 1))
acf(posterior$samples$chain2[, "lambda"]) # plot autocorrelation of lambda sample - chain 2
par(mfrow = c(1, 2))
plot(density(posterior$samples$chain1[ , "lambda"]), type = "l", xlab = expression(alpha),
     ylab = "Posterior density")
plot(density(posterior$samples$chain1[ , "averagelifetime"]), type = "l", xlab = expression(beta),
     ylab = "Posterior density")
posterior$summary$all.chains