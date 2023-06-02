## FUNCTIONS DEFINITIONS ##

################################################################################

lower_bounds = function(theta, x){
  
  lowb = matrix(rep(theta, dim(x)[1]), nrow=dim(x)[1], byrow=TRUE)
  lowb[x==0] = -Inf
  
  return(lowb)
}

################################################################################

upper_bounds = function(theta, x){
  
  lowb = matrix(rep(theta, dim(x)[1]), nrow=dim(x)[1], byrow=TRUE)
  lowb[x==1] = Inf
  
  return(lowb)
}

################################################################################

# Sigma --> covariance matrix of the Gaussian distribution;
# lower --> n*q matrix which stores in each row the lower bound of every z_i;
# upper --> n*q matrix which stores in each row the upper bound of every z_i;
# algorithm --> the algorithm we use to sample each z_i.

generate_Z = function(Sigma, lower, upper, algorithm="gibbs"){
  
  Z = matrix(NA, nrow=dim(lower)[1], ncol=dim(lower)[2])
  
  for(i in 1:dim(lower)[1]){
    Z[i,]=rtmvnorm(n=1, sigma=Sigma, lower=lower[i,], upper=upper[i,], algorithm=algorithm)
  }
  
  return(Z)
}

################################################################################

# Compute a step of Metropolis-Hastings algorithm to update the cutoffs vector.

MH_theta = function(Sigma, x, tau_prior, theta, sd_proposal=2){
  
  for(j in 1:dim(x)[2]){
    prop_theta_j = proposal_theta_j(theta[j], sd_proposal)
    temp_theta = theta
    temp_theta[j] = prop_theta_j
    
    acceptance = min(1, exp( log_density_theta(temp_theta, Sigma, x, tau_prior) - log_density_theta(theta, Sigma, x, tau_prior)))
    
    if(rbern(1,acceptance)){
      theta[j] = prop_theta_j
    }
  }
  
  return(theta)
}

################################################################################

# Sample from the proposal q(.|theta_j_old) ~ N(theta_j_old,sd)

proposal_theta_j = function(theta_j_old, sd=1){
  
  return(rnorm(n=1, mean=theta_j_old, sd=sd))
}

################################################################################

# Compute the density of the proposal q(theta_new_j|theta_j) ~ N(theta_j, sd^2)

dproposal = function(theta_new_j, theta_j, sd=1, log=FALSE){
  
  return(dnorm(x=theta_new_j, mean=theta_j, sd=sd, log=log))
}

################################################################################

# theta --> q dimensional vector for which we want to evaluate the density;
# Sigma --> q*q covariance matrix;
# x --> n*q matrix with 0 and 1 corresponding to the categorical observations;
# tau_prior --> prior variance of the cutoffs.

log_density_theta = function(theta, Sigma, x, tau_prior){
  
  lower_bounds_matrix = lower_bounds(theta, x)
  upper_bounds_matrix = upper_bounds(theta, x)
  
  log_dens = 0
  
  for(i in 1:dim(x)[1]){
    log_dens = log_dens + log(pmvnorm(lower=lower_bounds_matrix[i,], upper=upper_bounds_matrix[i,], sigma=Sigma)) 
  }
  
  for(j in 1:dim(x)[2]){
    log_dens = log_dens + dnorm(theta[j], sd=tau_prior, log=TRUE)
  }
  
  return(log_dens)
}

################################################################################
