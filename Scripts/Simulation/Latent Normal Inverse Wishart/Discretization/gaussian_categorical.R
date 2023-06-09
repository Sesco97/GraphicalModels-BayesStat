# Source useful scripts 

source("../../../Utilities/utilityFunctions.R")
source("../../../Utilities/cutoff_utility.R")
source("../Gaussian/gaussian.R")

# Load useful libraries

library(CholWishart)
library(tmvtnorm)

## FUNCTIONS DEFINITIONS ##

################################################################################

# Metropolis-Hastings for categorical Gaussian.

# We need to add latent Gaussian Z and the cutoffs vector theta.
# Theta will be sampled with a MH step.
# Z will be sampled from a truncated Gaussian according to the cutoffs theta.
# Once we have z we just follow the same procedure as for the pure Gaussian case.

MetropolisHastingsGaussianCategorical = function(data, initialCandidate, n.iter, burnin = 0, thin = 1, prior, p = NULL, b = NULL, algorithm="gibbs"){
  
  # We check that the passed parameters are correct
  
  if(!prior %in% c("Uniform","Binomial","Beta-Binomial")){
    stop("prior should be either 'Uniform', 'Binomial' or 'Beta-Binomial'!")
  }
  
  if(!isDecomposable(initialCandidate)){
    stop("Initial candidate graph should be decomposable!")
  }
  
  currentCandidate = initialCandidate
  
  b = 1 # degrees of freedom of the Hyper Inverse Wishart
  n = dim(data)[1]
  D = 10 * diag(1, dim(data)[2]) #parameter D of the Hyper Inverse Wishart
  x = data.matrix(data)
  tau_prior = 10 # std deviation prior on the cutoffs vector
  
  # Initialization
  
  theta = rep(0, dim(x)[2]) # initialize theta with zeros
  Z = generate_Z(Sigma=D, lower=lower_bounds(theta,x), upper=upper_bounds(theta,x), algorithm=algorithm)
  Sigma = D + t(Z)%*%Z
  D_star = Sigma # that's the parameter of the posterior for the Sigma
  
  # Run the burn-in iterations
  
  if(burnin!=0){
    message("BURN-IN")
    progressBarBI = txtProgressBar(min = 0, max = burnin, initial = 0, style = 3) 
    
    for(i in 1:burnin){
      setTxtProgressBar(progressBarBI,i)
      
      theta = MH_theta(Sigma=Sigma, x=x, tau_prior=tau_prior, theta=theta) # update theta vector with a MH step
      Sigma = update_Sigma(df=b+n, Dstar=D_star, adj=currentCandidate) # update Sigma
      Z = generate_Z(Sigma=Sigma, lower=lower_bounds(theta,x), upper=upper_bounds(theta,x), algorithm=algorithm)
      D_star = D + t(Z)%*%Z # update D_star for the HWS
      
      newCandidate = newGraphProposal(currentCandidate)
      num = logMarginalLikelihoodGaussian(newCandidate, Z, b, D, D_star=D_star)
      den = logMarginalLikelihoodGaussian(currentCandidate, Z, b, D, D_star=D_star)
      marginalRatio = exp(num - den)
      priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = betaBinomialPrior(currentCandidate,newCandidate,a,b))
      acceptanceProbability = min(marginalRatio * priorRatio,1)
      accepted = rbern(1,acceptanceProbability)
      
      if(accepted == 1){
        currentCandidate = newCandidate
      }
    }
    
    close(progressBarBI)
  }
  
  # Run the chain
  
  message("Metropolis-Hastings")
  progressBar = txtProgressBar(min = 0, max = n.iter, initial = 0, style = 3) 
  chain = list()
  c = 0
  
  for(i in 1:n.iter){
    setTxtProgressBar(progressBar,i)
    
    theta = MH_theta(Sigma=Sigma, x=x, tau_prior=tau_prior, theta=theta) # update theta vector with a MH step
    Sigma = update_Sigma(df=b+n, Dstar=D_star, adj=currentCandidate) # update Sigma
    Z = generate_Z(Sigma=Sigma, lower=lower_bounds(theta,x), upper=upper_bounds(theta,x), algorithm="gibbs")
    D_star = D + t(Z)%*%Z # update D_star for the HWS
    
    newCandidate = newGraphProposal(currentCandidate)
    num = logMarginalLikelihoodGaussian(newCandidate, Z, b, D, D_star=D_star)
    den = logMarginalLikelihoodGaussian(currentCandidate, Z, b, D, D_star=D_star)
    marginalRatio = exp(num - den)
    priorRatio = switch(prior, "Uniform" = 1, "Binomial" = binomialPrior(currentCandidate,newCandidate,p), "Beta-Binomial" = betaBinomialPrior(currentCandidate,newCandidate,a,b))
    acceptanceProbability = min(marginalRatio * priorRatio,1)
    accepted = rbern(1,acceptanceProbability)
    c = c + accepted
    
    if(accepted == 1){
      currentCandidate = newCandidate
    }
    
    if(i %% thin == 0){
      chain[[i/thin]] = currentCandidate
    }
  }
  
  if(burnin!=0){
    close(progressBarBI)
  }
  
  cat(paste0("\nThe average acceptance rate is: ", as.character(c / n.iter),"\n"))
  
  return(chain)
}

################################################################################

update_Sigma = function(df, Dstar, adj){
  
  inv.covariance = rgwish(n = 1, adj = adj, b = df, D = Dstar)
  covariance = solve(inv.covariance)
  
  return(covariance)
}

################################################################################

