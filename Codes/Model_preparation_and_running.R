#The model code below is written program R and uses the R2jags package to run JAGS.
#The data file, 'D', includes columns representing the study site, species (SppCode), number of 
#days during which the species was photographed ???1 time (detections), number of sampling 
#occasions (num.nights), covariate values, and diet classification (Diet)
#Define the group variables
D=read.csv("Input/MSOM_datafile.csv")%>%
  select(-X)
G <- cbind(as.numeric(D$Diet==1),as.numeric(D$Diet==0))
#Define the covariates for occupancy
X = select(D,sitecov1)
#Define the group covariates
XG = cbind(X*G[,1],X*G[,2])
#Define the covariates for detection
dX = select(D,obsvcov1:obsvcov2)
#Load the necessary libraries
library(R2jags); library(reshape2); library(dplyr)
#Define the necessary arguments to run the jags command
#Load all the data including the detection array, number of sampling occasions, individual 
#species sampled, total number of sampled species, and covariate information
data <- list(D = D$detection, N = ceiling(D[,"num.samples"]), Species =D$Sppcode, 
             n = nrow(D), nspp = max(as.numeric(D$Sppcode)),X = X, XG = XG,dX = dX)
#Specify the initial values
inits = function() {list(Z = as.numeric(data$D>0))}
#Specify the parameters to be monitored
params = c("rho","pbeta","spbeta","sigpbeta","mbeta","sigbeta","sbeta","gbeta",
           "psi.mean","sigma.occ","p.mean","sigma.p","alpha","Z","P")
#Specify the number of chains (nc), number of iterations (ni), burn-in period (nb), and thinning 
#rate (nthin)
nc = 2
ni = 6000
nb = 1000
nthin = 50#Write the model code to a text file called "AllMammals.txt"
cat(
  "    data {
      nBeta <- dim(X)
      nG <- dim(XG)
      nP <- dim(dX)
    }
    model {
      # Define covariance parameter between detection and mean occupancy
      rho ~ dunif(-1,1)
      var.p <- sigma.p /(1.-pow(rho,2))
      
      #Define prior distributions for occupancy parameters
      alpha.mean <- log(psi.mean) - log(1-psi.mean)
      psi.mean ~ dunif(0,1)
      sigma.occ ~ dunif(0,10)
      tau.occ <- pow(sigma.occ,-2)
      #Define prior distributions for true positive detections
      p.mean ~ dunif(0,1)
      b <- log(p.mean) - log(1-p.mean)
      sigma.p ~ dunif(0,10)
      tau.p <- pow(sigma.p,-2)
      #Define prior distributions for occupancy effects where nbeta is the number of occupancy 
      #covariates in the model, mbeta is the community-level hyper-parameter for each of the nbeta 
      #covariates, tbeta is the amount of variability in each of the community-level hyper-parameters, 
      #and sbeta is the species-specific covariate effects
      for (a in 1:nBeta[2]){
        mbeta[a] ~ dnorm(0,0.01) 
        sigbeta[a] ~ dunif(0,10)
        tbeta[a] <- pow(sigbeta[a],-2) 
        for (i in 1:(nspp)) { 
          sbeta[i,a] ~ dnorm(0,tbeta[a]) 
        }
      }
      #Define prior distributions for nG, the group-level hyper-parameters
      for (a in 1:nG[2]){
      gbeta[a] ~ dnorm(0,0.01)
    }
    #Define prior distributions for detection effects where nP is the number of detection covariates in 
    #the model, pbeta is the community-level hyper-parameter for each of the nP covariates, tpbeta is 
    #the amount of variability in each of the community-level hyper-parameters, and spbeta is the 
    #species-specific covariate effects
    for (a in 1:nP[2]){
      pbeta[a] ~ dnorm(0,0.01)
      sigpbeta[a] ~ dunif(0,10)
      tpbeta[a] <- pow(sigpbeta[a],-2) 
      for (i in 1:(nspp)) { 
        spbeta[i,a] ~ dnorm(0,tpbeta[a]) 
      }
    }
    #Define prior distributions for the occupancy and detection covariates for each species 
    for (i in 1:(nspp)) {
      alpha[i] ~ dnorm(alpha.mean, tau.occ)
      mu.p[i] <- b + (rho*sigma.p /sigma.occ)*(alpha[i] - alpha.mean)
      P[i] ~ dnorm(mu.p[i], var.p)
    }
    
    #Estimate the occupancy probability (latent Z matrix) for each species at each camera station
    for (j in 1:n) {
      logit(psi[j]) <- alpha[Species[j]] + inprod(mbeta,X[j,])+ inprod(sbeta[Species[j],],
                                                                       X[j,]) + inprod(gbeta,XG[j,])
      #Estimate the detection probability for each species at each camera station
      logit(p[j]) <- P[Species[j]] + inprod(pbeta,dX[j,]) + inprod(spbeta[Species[j],],
                                                                   dX[j,]) 
      Z[j] ~ dbern(psi[j])
      zp[j] <- p[j]*Z[j]
      D[j] ~ dbin(zp[j], N[j])
    }
    }
", file = "Input/JAGS_model.txt")

#Run the model and call the results "output"
output <- jags(data = data, inits = inits, parameters.to.save = params, 
               model.file ="Input/JAGS_model.txt", 
               n.chains =nc, n.iter =ni, n.burnin =nb, n.thin =nthin)
write_rds(output,"Results/model_output.Rdata")

