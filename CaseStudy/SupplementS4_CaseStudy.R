##
##
## Case study described in Hostetter NJ, NJ Lunn, ES Richardson, EV Regehr, and SJ Converse: 
## Age-structured Jolly-Seber model expands inference and improves parameter estimation from capture-recapture data.
##
## This script loads, formats, and runs the models described in the case study. 
## For detailed examples of data and models, we recommend the simulation script (SupplementS1_SimulationAgeStructuredJS.R)
##
## models include
## JS without age data ("NoAge")
## JS with age data and constant survival ("Age")
## JS with age data and survival is a quadratic function of age ("AgeSurvival")
##
## 10-July-2020
##
## Although these data have been processed successfully on a computer system at the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the display or utility of the data for other purposes, nor on all computer systems, nor shall the act of distribution constitute any such warranty. The USGS or the U.S. Government shall not be held liable for improper or incorrect use of the data described and/or contained herein.

## required packages
library(nimble)


## load data
dat <- dget("SupplementS3_CaseStudyData.txt")
list2env(dat, .GlobalEnv) # send list variable to global environment for ease of indexing


# ydatAug: capture-recapture data, inlcuding augmented individuals 
# n: number of detected individuals
# nz: number of augmented individuals (note M = n+nz)
# M: augmentation value
# K: study occasions
# ageMatAugMinusOne: age matrix, including augmented individuals. Ages are 1-actual ages as individuals *enter* at actual age 2 (independence)
# centered_age: median observed age (minus 1)
# max.age: maximum age in occasion 1 (minus 1).
# zdatNoAGE: known values of z in JS model without age data
# zstNoAGE: initial values for z in JS model without age data. Prevents initialization errors
# zdatAGE: known values of z in JS model with age data
# zstAGE: initial values for z in JS model with age data. Prevents initialization errors


## MCMC settings (warning, hours of run time with full settings)
nc <- 6; nAdapt <- 500; nb <- 20000; ni <- 200000+nb; nthin <- 10
nc <- 2; nAdapt <- 500; nb <- 100; ni <- 500+nb; nthin <- 1 # demonstrate model

## select model to run
run.model <- "NoAge" # "NoAge", "AgeSurvival","Age" 



if(run.model=="NoAge"){

##
##
##
## OPTION 1: no age data
##
##
##
##
  
  code_JS <- nimbleCode({
    
    # Priors and constraints
    psi ~ dbeta(1,1)          # data augmentation
    p ~ dbeta(1,1)            # detection probability
    phi ~ dbeta(1,1)          # survival
 
    # recruitment prob at k 
    beta[1:K] ~ ddirch(b[1:K])

    eta[1] <- beta[1]
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }


    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)

      # At first occasion
      # State process
      u[i,1] ~ dbern(eta[1])   
      z[i,1] <- u[i,1]*w[i]

      # Observation process
      y[i,1] ~ dbern(z[i,1]*p)
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for 'recruitment'. 
      recruit[i,1] <- z[i,1]             # 'recruited' at k
 
      # for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi + avail[i,t-1]*eta[t])   
        z[i,t] <-u[i,t]*w[i]   

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])            # still available for 'recruitment'. 
        recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # 'recruited' at k
      } #t
    } #i
    
    #derived population level stuff
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               # Annual abundance
      B[t] <- sum(recruit[1:M,t])         # Number of entries
    } #t

    Nsuper <- sum(w[1:M])            # Superpopulation size

    # Optional Freeman-Tukey GOF evaluating observed and expected counts.
    for(t in 1:K){
     expectedCount[t] <- N[t]*p
     Count.new[t] ~ dbin(p, N[t])

     ## discrepancy from original
     E.dat[t] <- (sqrt(Counts[t])-sqrt(expectedCount[t]) )^2

     ## discrepancy from new
     E.new[t] <- (sqrt(Count.new[t])-sqrt(expectedCount[t]) )^2
    }
    fit.dat <- sum(E.dat[1:K])
    fit.new <- sum(E.new[1:K])
    
  })#END model




## Format for NIMBLE without age
nim.data<- list(u=zdatNoAGE,y = ydatAug, w=c(rep(1, n), rep(NA, nz)), b=rep(1,K)  ) 
nim.constants<- list(K=K, M=M, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi", "beta", "psi", "Nsuper", "N","fit.dat", "fit.new")
inits <- list(phi = .80, p = .2, u = zstNoAGE , w=c(rep(NA, n), rep(0, nz)), psi = .5, beta=rep(1/K, K) )


## RUN!
  ## Nimble steps
  Rmodel <- nimbleModel(code=code_JS, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
  conf <- configureMCMC(Rmodel,monitors=params, control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
  Rmcmc <- buildMCMC(conf)  #produces an uncompiled R mcmc function
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
  
}#if(run.model=="NoAge")




if(run.model=="Age"){
##
##
##
## OPTION 2: with age data, constant survival
##
##
##
##
  
  ## NIMBLE MODEL
  
  code_JS_AGE <- nimbleCode({
    
    # Priors and constraints
    psi~ dbeta(1,1)           # data augmentation
    p ~ dbeta(1,1)            # detection
    phi ~ dbeta(1,1)          # survival
   
    # recruitment prob at k 
    beta[1:K] ~ ddirch(b[1:K])

    eta[1] <- beta[1]
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }

    # starting age distribution is not yet recruited (age = 0), or age if recruited
    piAGE[1:max.age] ~ ddirch(a[1:max.age])
    piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  
 
     
    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)

      # Define latent state at first occasion
      # Age process
      agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
      age[i,1] <- agePlusOne[i]-1
      trueAge[i,1] <- age[i,1]*z[i,1]    # age if alive (zero if not yet entered or dead)

      # State process
      u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] >0
      z[i,1] <- u[i,1]*w[i]  

      # Observation process
      y[i,1] ~ dbern(z[i,1]*p)
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for recruitment. 
      recruit[i,1] <- z[i,1]             # recruited at k
 
      # for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi + avail[i,t-1]*eta[t])   
        z[i,t] <- u[i,t]*w[i]  
        trueAge[i,t] <- age[i,t]*z[i,t]    
       
        # Age process
        age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])       
        recruit[i,t] <- equals(z[i,t]-z[i,t-1],1) # recruited at k
      } #t
    } #i
    
    #derived population level stuff
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               # Annual abundance
      B[t] <- sum(recruit[1:M,t])         # Number of entries
    } #t

    Nsuper <- sum(w[1:M])            # Superpopulation size

    # Optional Freeman-Tukey GOF evaluating observed and expected counts.
    for(t in 1:K){
     expectedCount[t] <- N[t]*p
     Count.new[t] ~ dbin(p, N[t])

     ## discrepancy from original
     E.dat[t] <- (sqrt(Counts[t])-sqrt(expectedCount[t]) )^2

     ## discrepancy from new
     E.new[t] <- (sqrt(Count.new[t])-sqrt(expectedCount[t]) )^2
    }
    fit.dat <- sum(E.dat[1:K])
    fit.new <- sum(E.new[1:K])
    
  })#END model
  

##
## Format for NIMBLE with age
nim.data<- list(u=zdatAGE,y = ydatAug, age=ageMatAugMinusOne, agePlusOne=ageMatAugMinusOne[,1]+1, w=c(rep(1, n), rep(NA, nz)),
                b=rep(1,K), a=rep(1, max.age) ) 
nim.constants<- list(K=K, M=M, max.age=max.age, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi", "beta", "psi", "piAGE", "Nsuper", "N","fit.dat", "fit.new")   
inits <- list(phi = .9, p = .2, psi=.6, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), 
             beta=rep(1/K, K),agePlusOne=c(rep(NA, n), rep(1, nz)) )  


## RUN!
  ## Nimble steps
  Rmodel <- nimbleModel(code=code_JS_AGE, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
  conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
  Rmcmc <- buildMCMC(conf)  
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)
 
}#if(run.model=="Age")





if(run.model=="AgeSurvival"){
##
##
##
##
## Option 3. with age data, survival a quadratic function of age 
##
##
##
##




## NIMBLE MODEL

code_JS_AgeSurvival <- nimbleCode({

# Priors and constraints
    psi~ dbeta(1,1)              # data augmentation 
    p ~ dbeta(1,1)               # detection
    phi0 ~ dbeta(1,1)            # survival at centered age
    alpha0 <- logit(phi0)
    alpha1 ~ dnorm(0, sd=10)     # relationship with age[i,t]
    alpha2 ~ dnorm(0, sd=10)     # relationship with age[i,t]^2

    # recruitment prob at k 
    beta[1:K] ~ ddirch(b[1:K])

    eta[1] <- beta[1]
    for(k in 2:K){
     eta[k] <- beta[k]/(1-sum(beta[1:(k-1)]))
    }

    # starting age distribution is not yet recruited (age = 0), or age if alive
    piAGE[1:max.age] ~ ddirch(a[1:max.age])
    piAGEuncond[1:(max.age+1)] <- c( (1-eta[1] ), eta[1]*piAGE[1:max.age] )  

    # Likelihoods 
    for (i in 1:M){
      w[i] ~ dbern(psi)

      # Define latent state at first occasion
      # Age process
      agePlusOne[i] ~ dcat(piAGEuncond[1:(max.age+1)]) 
      age[i,1] <- agePlusOne[i]-1
      trueAge[i,1] <- age[i,1] * z[i,1]
  
      # State process
      u[i,1] <- step(agePlusOne[i]-1.1)  # alive if age[i,1] >0
      z[i,1] <- u[i,1]*w[i]  

      # Observation process
      y[i,1] ~ dbern(z[i,1]*p)
      
      # derived stuff
      avail[i,1] <- 1- u[i,1]            # still available for recruitment
      recruit[i,1] <- z[i,1]             # recruited at t
 
      # for occasions >1     
      for (t in 2:K){

        # State process
        u[i,t] ~ dbern(u[i,t-1]*phi[i,t-1] + avail[i,t-1]*eta[t])   
        logit(phi[i,t-1]) <- alpha0 + alpha1*(age[i,t-1]-centered_age) + alpha2*(age[i,t-1]-centered_age)^2        # centered
        z[i,t] <- u[i,t]*w[i]  

        # Age process
        age[i,t] <- age[i,t-1] + max(u[i,1:t]) # ages by one year after recruitment 
        trueAge[i,t] <- age[i,t] * z[i,t]

        # Observation process
        y[i,t] ~ dbern(z[i,t]*p)

        # derived stuff
        avail[i,t] <- 1- max(u[i,1:t])       
        recruit[i,t] <- equals(z[i,t]-z[i,t-1],1)  # recruited at t
      } #t
    } #i
    
    #derived population level stuff
    for (t in 1:K){
      N[t] <- sum(z[1:M,t])               # Annual abundance
      B[t] <- sum(recruit[1:M,t])         # Number of entries
    } #t

    Nsuper <- sum(w[1:M])            # Superpopulation size

   # monitor age specific survival (buffer beyond expected max.age if interested). This could be done outside the model.   
   for(aa in 1:(max.age+5)){ 
    logit(phi.age[aa]) <- alpha0 + alpha1*(aa-centered_age) + alpha2*(aa-centered_age)^2 
   }    

    # Optional Freeman-Tukey GOF evaluating observed and expected counts.
    for(t in 1:K){
     expectedCount[t] <- N[t]*p
     Count.new[t] ~ dbin(p, N[t])

     ## discrepancy from original
     E.dat[t] <- (sqrt(Counts[t])-sqrt(expectedCount[t]) )^2

     ## discrepancy from new
     E.new[t] <- (sqrt(Count.new[t])-sqrt(expectedCount[t]) )^2
    }
    fit.dat <- sum(E.dat[1:K])
    fit.new <- sum(E.new[1:K])

  })#END model


##
## Format for NIMBLE 

nim.data<- list(u=zdatAGE,y = ydatAug, age=ageMatAugMinusOne, agePlusOne=ageMatAugMinusOne[,1]+1, w=c(rep(1, n), rep(NA, nz)),
                b=rep(1,K), a=rep(1, max.age)) 
nim.constants<- list(K=K, M=M, max.age=max.age, centered_age=centered_age, Counts=apply(ydatAug, 2,sum))
params <- c("p", "phi.age", "phi0","alpha0", "alpha1", "alpha2",
              "beta", "piAGE","psi", "Nsuper", "N","fit.dat", "fit.new") # save trueAge to derive estimated age structure our side the model.
inits <- list(phi0 = .75, alpha1 = 0, alpha2 = 0, p = .2, psi=.6, u = zstAGE, w=c(rep(NA, n), rep(0, nz)), 
              beta=rep(1/K, K),agePlusOne=c(rep(NA, n), rep(1, nz)) )


## RUN!
  ## Nimble steps
  Rmodel <- nimbleModel(code=code_JS_AgeSurvival, constants=nim.constants, data=nim.data, calculate=F, check=F, inits=inits)
  conf <- configureMCMC(Rmodel,monitors=params,control = list(adaptInterval = nAdapt), thin=nthin, useConjugacy = TRUE)
  Rmcmc <- buildMCMC(conf)  
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  out  <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits=inits, thin=nthin,
                 setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)

}# if(run.model=="AgeSurvival")



summary(out)
