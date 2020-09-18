# Alternative coding of the MCMC sampler.

# Not sure what this function is (possibly the likelihood for psi?)
target_yq <- function(psi,weights,MBM,tau){
  sum(log(weights))-0.5*t(psi)%*%MBM%*%psi*tau
}

target_beta <- function(beta,weights,beta_prior_mean,g,tau){
  sum(log(weights))-0.5*(t(beta-beta_prior_mean)%*%(beta-beta_prior_mean))*g*tau
}

yq_giv_tau <- function(psi,tau,MBM){
  0.5 * length(psi) * log(tau) -0.5*(t(psi)%*%MBM%*%psi)*tau
  }

beta_giv_tau<-function(beta,beta_prior_mean,tau,g){
  btb<-as.vector(t(beta-beta_prior_mean)%*%(beta-beta_prior_mean))
  btb_tau<-btb*g*tau
  log(g*tau)*length(beta)/2-0.5*btb_tau
}

tau_prior <- function(tau,tau_alpha_1,tau_alpha_2){(1+tau_alpha_1)*log(tau)-(tau_alpha_2/tau)}

target_tau<-function(tau,beta,beta_prior_mean,tau_alpha_1,tau_alpha_2,psi,MBM,g){
  p <- length(beta)
  q <- nrow(MBM)
  # ((p+q)/2)*log(tau)+ 
    yq_giv_tau(psi,tau,MBM)+ 
    beta_giv_tau(beta,beta_prior_mean,tau,g)+
    tau_prior(tau,tau_alpha_1,tau_alpha_2)
}

move_tau <- function(tau, tau_proposed,beta_prior_mean,tau_alpha_1,tau_alpha_2,psi,MBM,g){
  # Compute the probability density ratio:
  pdr_tau <- target_tau(tau_proposed,beta,beta_prior_mean = beta_prior_mean,tau_alpha_1,tau_alpha_2,psi,MBM,g = g)-
    target_tau(tau,beta,beta_prior_mean = beta_prior_mean,tau_alpha_1,tau_alpha_2,psi,MBM,g = g)
  # message("Tau PDR: ",round(exp(pdr_tau),3)) ##DEBUG##
  # Accept or reject proposed value of beta
  rexp(1) > -(pdr_tau + log(tau_proposed)-log(tau))
}



propose_beta <- function(beta,cov_mat){
  mvrnorm(1,beta,cov_mat)
}

propose_tau <- function(tau,sd){
  exp(rnorm(1,mean= log(tau),sd=sd))
}

propose_psi <- function(psi,sd){
  z <- rnorm(length(psi),0,sd)
  c(psi + inv_chol_MBM %*% z)
}

move_beta <- function(beta,beta_proposed,beta_prior_mean,weights,weights_proposed,g = 10,tau){
  # Compute the probability density ratio:
    pdr_beta <- target_beta(beta = beta_proposed,weights=weights_proposed,beta_prior_mean = beta_prior_mean,g=g,tau)-
    target_beta(beta = beta,weights=weights,beta_prior_mean = beta_prior_mean,g=g,tau)
    # message("Beta PDR: ",round(exp(pdr_beta),3)) ##DEBUG##
  # Accept or reject proposed value of beta
    rexp(1) > -pdr_beta
}

move_psi <- function(psi,psi_proposed,weights,weights_proposed,MBM, tau){
  # Compute the probability density ratio:
  pdr_psi <- target_yq(weights=weights_proposed, psi=psi_proposed, MBM, tau=tau)-
               target_yq(w=weights, psi=psi, MBM, tau=tau)
  # message("Psi PDR: ",round(exp(pdr_psi),3)) ##DEBUG##
  # Accept or reject proposed value of psi
  rexp(1) > -pdr_psi
}

# Simulate some data:
# {
#   tau <- 0.2
#   psi <- mvrnorm(n = 1, mu = rep(0,nrow(MBM)), Sigma = tau * MBM)
#   beta <- c(0,0.1)
#   y <- x %*% beta + M %*% psi
# }

# plot(beta_sample[2,], apply(psi_sample,2,mean))

niter = 10000
sd_beta <- 0.01
sd_tau <- 0.8
sd_psi <- 0.4
g = 10
beta_prior_mean <- c(0,0)
tau_alpha_1 <- 1
tau_alpha_2 <- 1

# Compute values needed:
chol_MBM <- chol(MBM)
inv_chol_MBM <- solve(chol_MBM)
# Initialise parameter values:
tau = 10
psi = rnorm(nrow(MBM))
# beta = rep(0,ncol(x))
beta = c(0,0.1)
# Create objects to store results:
psi_sample <- matrix(0,nrow= nrow(MBM), ncol=niter)
beta_sample <- matrix(0,nrow= ncol(x), ncol=niter)
tau_sample <- rep(NA,niter)
el_sample <- rep(NA, niter)
# Initialise weights for current parameter values:
mu <-x %*% beta + M %*% psi
weights<-el.test(y-mu, 0)$wts
weights<-weights/sum(weights)

for(iter in 1:niter){
  
  # message("Iteration: ",iter) ##DEBUG##
  
  #### PSI ######################

  # Propose a new value of psi:
  psi_proposed <- propose_psi(psi,sd = sd_psi / sqrt(tau))
  # Compute new weights:
  mu_proposed  <-x %*% beta + M %*% psi_proposed
  weights_proposed <- el.test(y-mu_proposed, 0)$wts
  weights_proposed <- weights_proposed/sum(weights_proposed)
  # Accept or Reject proposed psi:
  MOVE_PSI <- move_psi(psi,psi_proposed,weights,weights_proposed,MBM, tau)
  if(MOVE_PSI){
    psi <- psi_proposed
    weights <- weights_proposed
    mu <- mu_proposed
    # message("Psi accepted.") ##DEBUG##
  }else{
    # message("Psi rejected.") ##DEBUG##
  }
  
  # #### BETA #####################
  # Propose a new value of beta:
  beta_proposed <- propose_beta(beta,cov_mat = sd_beta*diag(length(beta)))
  # Compute new weights:
  mu_proposed  <-x %*% beta_proposed + M %*% psi
  weights_proposed <- el.test(y-mu_proposed, 0)$wts
  weights_proposed <- weights_proposed/sum(weights_proposed)
  # Accept or Reject proposed psi:
  MOVE_BETA <- move_beta(beta,beta_proposed,beta_prior_mean,weights,weights_proposed,g = g,tau)
  if(MOVE_BETA){
    beta <- beta_proposed
    weights <- weights_proposed
    mu <- mu_proposed
    # message("Beta accepted.") ##DEBUG##
  }else{
    # message("Beta rejected.") ##DEBUG##
  }
  
  #### TAU ######################
  # Propose a new value of TAU:
  tau_proposed <- propose_tau(tau,sd_tau)
  # Compute new weights:
  MOVE_TAU <- move_tau(tau, tau_proposed,beta_prior_mean,tau_alpha_1,tau_alpha_2,psi,MBM,g)
  if(MOVE_TAU){
    tau <- tau_proposed
    # message("Tau accepted.") ##DEBUG##
  }else{
    # message("Tau rejected.") ##DEBUG##
  }
  
  #### SAVE VALUES ##############
  psi_sample[,iter] <- psi
  beta_sample[,iter] <- beta
  tau_sample[iter] <- tau
  el_sample[iter] <- sum(log(weights))
}



par(mfrow = c(3,2))
plot(beta_sample[1,],type = "l")
plot(beta_sample[2,],type = "l")
plot(tau_sample,type = "l")
plot(psi_sample[2,],type = "l")
plot(el_sample,type = "l")
plot(y~ I(x %*% apply(beta_sample,1,mean) + M %*% apply(psi_sample,1,mean)))
abline(0,1)

acceptance_rate <- function(x){
  1-mean(diff(x) == 0)
}



lm(y ~ x[,2])
apply(beta_sample,1,mean)
# plot(y~mu_proposed)
# abline(0,1)

acceptance_rate(beta_sample[1,])
acceptance_rate(psi_sample[1,])
acceptance_rate(tau_sample)
# 

