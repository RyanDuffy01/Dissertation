#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(nimble)
library(igraph)
library(coda)
library(R6)

# Rats example. Random effects model 

################################################################
################################################################

# Specify the statistical model
rats_reCode <- nimbleCode({
  
  # Specify the likelihood:
  
  for(i in 1 : N) {
    for(j in 1 : T) {
      Y[i,j] ~ dnorm(mu[i,j],tau)
      mu[i,j] <- alpha[i] + beta*(x[j] - mean(x[]))/sd(x[])
    }		
    alpha[i] <- alpha_mean + epsilon[i]
  }
  
  # Prior specification:
  
  alpha_mean ~ dnorm(0,0.00001)
  beta ~ dnorm(0,0.00001)
  tau ~ dgamma(0.001,0.001)
  sigma2 <- 1 / tau
  for(k in 1:N){
    epsilon[k] ~ dnorm (0, invkappa)
  }
  kappa<- 1/invkappa
  invkappa ~ dgamma(0.001, 0.001)
  
})

# Values for some constants in the model
rats_reConsts <- list(N = 30, T=5)

# The data values
rats_reData <- list(x = c(8.0, 15.0, 22.0, 29.0, 36.0),	
                  Y = t(matrix( c(151, 199, 246, 283, 320, #NOTE THE DIFFERENCE TO OpenBUGS for
                                  145, 199, 249, 293, 354,   # reading data in matrix form!
                                  147, 214, 263, 312, 328,
                                  155, 200, 237, 272, 297,
                                  135, 188, 230, 280, 323,
                                  159, 210, 252, 298, 331,
                                  141, 189, 231, 275, 305,
                                  159, 201, 248, 297, 338,
                                  177, 236, 285, 350, 376,
                                  134, 182, 220, 260, 296,
                                  160, 208, 261, 313, 352,
                                  143, 188, 220, 273, 314,
                                  154, 200, 244, 289, 325,
                                  171, 221, 270, 326, 358,
                                  163, 216, 242, 281, 312,
                                  160, 207, 248, 288, 324,
                                  142, 187, 234, 280, 316,
                                  156, 203, 243, 283, 317,
                                  157, 212, 259, 307, 336,
                                  152, 203, 246, 286, 321,
                                  154, 205, 253, 298, 334,
                                  139, 190, 225, 267, 302,
                                  146, 191, 229, 272, 302,
                                  157, 211, 250, 285, 323,
                                  132, 185, 237, 286, 331,
                                  160, 207, 257, 303, 345,
                                  169, 216, 261, 295, 333,
                                  157, 205, 248, 289, 316,
                                  137, 180, 219, 258, 291,
                                  153, 200, 244, 286, 324),
                                5,30)))

# one set of initial values before building the model                 
rats_reInits <- list(alpha_mean=0,beta = 0, tau = 1, invkappa = 1) # Nimble will generate the rest 

# to build the model
rats_re <- nimbleModel(code = rats_reCode, name = "rats_re", constants = rats_reConsts,
                    data = rats_reData, inits<-rats_reInits)

# To compile the model
Crats_re <- compileNimble(rats_re)

# set up the monitored quantities. Default is all of the random quantities
rats_reConf <- configureMCMC(rats_re, monitors = c('alpha_mean','beta','kappa','sigma2','alpha'), print = TRUE) 
# build the MCMC algorithm
rats_reMCMC <- buildMCMC(rats_reConf)
# compile the MCMC chain 
Crats_reMCMC <- compileNimble(rats_reMCMC, project = rats_re)

####################################################################################
####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS  #
####################################################################################
set.seed(10)
rats_reInits <- list(list(alpha_mean=0,beta = 0, tau = 1, invkappa = 1), 
                     list(alpha_mean=1,beta = 1, tau = 2, invkappa = 2))
posterior <- runMCMC(Crats_reMCMC, niter = 5000, thin=1, nburnin=2000, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = rats_reInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
plot(combinedchains) # too many plots sometimes
plot(combinedchains[,c('alpha_mean','beta','kappa','sigma2')]) 
plot(combinedchains[,'alpha[1]'])
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
posterior$summary$all.chains

##################################################################################
##################################################################################
######### CODE TO PREDICT THE WEIGHT OF A RAT  (SEE LECTURE NOTES, SECTION 2.5)  #
##################################################################################

# Specify the statistical model
rats_predCode <- nimbleCode({
  
  # Specify the likelihood:
  
  for(i in 1 : N) {
    for(j in 1 : T) {
      Y[i,j] ~ dnorm(mu[i,j],tau)
      mu[i,j] <- alpha[i] + beta*(x[j] - mean(x[]))/sd(x[])
    }		
    alpha[i] <- alpha_mean + epsilon[i]
  }
  
  # Prior specification:
  
  alpha_mean ~ dnorm(0,0.00001)
  beta ~ dnorm(0,0.00001)
  tau ~ dgamma(0.001,0.001)
  sigma2 <- 1 / tau
  for(k in 1:N){
    epsilon[k] ~ dnorm (0, invkappa)
  }
  kappa<- 1/invkappa
  invkappa ~ dgamma(0.001, 0.001)
  
})

# Values for some constants in the model
rats_predConsts <- list(N = 31, T=5)

# The data values
rats_predData <- list(x = c(8.0, 15.0, 22.0, 29.0, 36.0),	
                    Y = t(matrix( c(151, 199, 246, 283, 320, #NOTE THE DIFFERENCE TO OpenBUGS for
                                    145, 199, 249, 293, 354,   # reading data in matrix form!
                                    147, 214, 263, 312, 328,
                                    155, 200, 237, 272, 297,
                                    135, 188, 230, 280, 323,
                                    159, 210, 252, 298, 331,
                                    141, 189, 231, 275, 305,
                                    159, 201, 248, 297, 338,
                                    177, 236, 285, 350, 376,
                                    134, 182, 220, 260, 296,
                                    160, 208, 261, 313, 352,
                                    143, 188, 220, 273, 314,
                                    154, 200, 244, 289, 325,
                                    171, 221, 270, 326, 358,
                                    163, 216, 242, 281, 312,
                                    160, 207, 248, 288, 324,
                                    142, 187, 234, 280, 316,
                                    156, 203, 243, 283, 317,
                                    157, 212, 259, 307, 336,
                                    152, 203, 246, 286, 321,
                                    154, 205, 253, 298, 334,
                                    139, 190, 225, 267, 302,
                                    146, 191, 229, 272, 302,
                                    157, 211, 250, 285, 323,
                                    132, 185, 237, 286, 331,
                                    160, 207, 257, 303, 345,
                                    169, 216, 261, 295, 333,
                                    157, 205, 248, 289, 316,
                                    137, 180, 219, 258, 291,
                                    153, 200, 244, 286, 324,
                                    150,NA,NA,NA,NA),
                                  5,31)))

# one set of initial values before building the model                 
rats_predInits <- list(alpha_mean=0,beta = 0, tau = 1, invkappa = 1) # Nimble will generate the rest 

# to build the model
rats_pred <- nimbleModel(code = rats_predCode, name = "rats_pred", constants = rats_predConsts,
                       data = rats_predData, inits<-rats_predInits)

# To compile the model
Crats_pred <- compileNimble(rats_pred)

# set up the monitored quantities. Default is all of the random quantities
rats_predConf <- configureMCMC(rats_pred, monitors = c('alpha_mean','beta','kappa','sigma2','Y'), print = TRUE) 
# build the MCMC algorithm
rats_predMCMC <- buildMCMC(rats_predConf)
# compile the MCMC chain 
Crats_predMCMC <- compileNimble(rats_predMCMC, project = rats_pred)

####################################################################################
####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS  #
####################################################################################
set.seed(10)
rats_predInits <- list(list(alpha_mean=0,beta = 0, tau = 1, invkappa = 1), 
                     list(alpha_mean=1,beta = 1, tau = 2, invkappa = 2))
posterior <- runMCMC(Crats_predMCMC, niter = 10000, thin=1, nburnin=5000, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = rats_predInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
plot(combinedchains[,c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]')]) #need the space after the comma
autocorr.plot(posterior$samples$chain1[,c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]')])
autocorr.plot(posterior$samples$chain2[,c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]')])
gelman.diag(combinedchains[,c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]')])
gelman.plot(combinedchains[,c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]')])
posterior$summary$all.chains[c('Y[31, 2]','Y[31, 3]','Y[31, 4]','Y[31, 5]'),]

# Note, you could unify the two sets of code above into one concise model
