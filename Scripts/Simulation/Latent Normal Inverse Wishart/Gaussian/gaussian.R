# Source useful scripts 

source("../../../Utilities/utilityFunctions.R")
source("../../../Utilities/cutoff_utility.R")

# Load useful libraries

library(CholWishart)
library(tmvtnorm)

## FUNCTIONS DEFINITIONS ##

################################################################################

# Generate a dataset of n.obs observations issued from a Multivariate Gaussian with dimension n.variables.
# Data is generated according to the graph given by adjacencyMatrix or by p, which is the vector of probability that an edge will be in the graph.
# If the covariance Sigma is null, the function will sample one from the Hyper Inverse Wishart distribution.

# The output is a list with the adjacency matrix, the data and the covariance which generated the data.

generateGaussianDataFromGraph = function(adjacencyMatrix=NULL, n.obs, n.variables, p=NULL, covariance=NULL){
  
  if(is.null(adjacencyMatrix)){
    
    while(TRUE){
      graph = erdos.renyi.game(n.variables,p,type="gnp",directed = FALSE)
      adjacencyMatrix = as_adjacency_matrix(graph, sparse = 0)
      
      if(isDecomposable(adjacencyMatrix)){
        break
      }
    }
  }
  
  if(!isDecomposable(adjacencyMatrix)){
    stop("Graph should be decomposable.")
  }
  
  if(is.null(covariance)){
    inv.covariance = rgwish(1, adj = adjacencyMatrix, D = 10 * diag(1,n.variables))
    covariance = solve(inv.covariance)
  }
  
  mu = c(rep(0, n.variables))
  data = dataCopy = data.frame(rmvnorm(n.obs, mu, covariance))
  
  return (list(adjacencyMatrix, data, covariance))
}

################################################################################

# Compute the log marginal likelihood for data and for the graph represented by adjacenyMatrix.
# If D_star is passed the code will be more efficient.

logMarginalLikelihoodGaussian = function(adjacencyMatrix, data, b, D, D_star=NULL){
  
  n = dim(data)[1]
  q = dim(data)[2]
  b_star = b + n
  
  if(is.null(D_star)){
    D_star = D + t(data.matrix(data))%*%data.matrix(data)
  }
  
  return (-log(2*pi)*n*q/2 + logh(adjacencyMatrix, b, D) - logh(adjacencyMatrix, b_star, D_star))
}

################################################################################

# Compute log(h) for the graph represented by adjacencyMatrix and parameter b and D.

logh = function(adjacencyMatrix, b, D){
  
  decomposition = getCliquesAndSeparators(adjacencyMatrix)
  cliques = decomposition[[1]]
  separators = decomposition[[2]]
  
  # Loop for each clique 
  
  num = 0
  
  for(i in 1:length(cliques)){
    clique = cliques[[i]]
    Dc = D[clique,clique]
    
    if(!is.matrix(Dc)){
      Dc = matrix(Dc)
    }
    
    cardC = length(clique)
    num = num + ((b + cardC - 1)/2)*determinant(Dc/2, logarithm=TRUE)$modulus[1] - lmvgamma((b + cardC - 1)/2, cardC) 
  }
  
  # Loop for each separator
  
  den = 0
  
  if(length(separators) == 0){
    den = 0
  }
  
  else{
    
    for(i in 1:length(separators)){
      separator = separators[[i]]
      Ds = D[separator, separator]
      
      if(!is.matrix(Ds)){
        Ds = matrix(Ds)
      }
      
      cardS = length(separator)
      den = den + ((b + cardS - 1)/2)*determinant(Ds/2, logarithm=TRUE)$modulus[1] - lmvgamma((b + cardS - 1)/2, cardS)
    }
  }
  
  result = num - den
  
  return(result)
}

################################################################################

MetropolisHastingsGaussian = function(data, initialCandidate, n.iter, burnin = 0, thin = 1, prior, p = NULL, b = NULL){
  
  # We check that the passed parameters are correct
  
  if(!prior %in% c("Uniform","Binomial","Beta-Binomial")){
    stop("prior should be either 'Uniform', 'Binomial' or 'Beta-Binomial'!")
  }
  
  if(!isDecomposable(initialCandidate)){
    stop("Initial candidate graph should be decomposable!")
  }
  
  currentCandidate = initialCandidate
  print(dim(data))
  b = 1 # degrees of freedom of the Hyper Inverse Wishart
  D = 10 * diag(1, dim(data)[2]) # parameter D of the Hyper Inverse Wishart
  D_star = D + t(data.matrix(data))%*%data.matrix(data)
  
  # Run the burn-in iterations
  
  if(burnin!=0){
    message("BURN-IN")
    progressBarBI = txtProgressBar(min = 0, max = burnin, initial = 0, style = 3) 
    
    for(i in 1:burnin){
      setTxtProgressBar(progressBarBI,i)
      newCandidate = newGraphProposal(currentCandidate)
      num = logMarginalLikelihoodGaussian(newCandidate,data, b, D, D_star=D_star)
      den = logMarginalLikelihoodGaussian(currentCandidate,data, b, D, D_star=D_star)
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
    newCandidate = newGraphProposal(currentCandidate)
    num = logMarginalLikelihoodGaussian(newCandidate,data,b, D, D_star=D_star)
    den = logMarginalLikelihoodGaussian(currentCandidate,data,b, D, D_star=D_star)
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
