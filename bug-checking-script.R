# creating a autoregression matrix randomly to model BYM and leroux 

library(MASS)
library(gtools)
# creating autoregression matrix randomly for a set of 50 small areas

R<- matrix(sample(-1:0,2500,replace=T),nrow=50, ncol=50)
# R<- matrix(0,nrow=50, ncol=50)
# R[cbind(2:50,1:49)] <- -1
# R[cbind(1:49,2:50)] <- -1
diag_R<-c()
for(i in 1:50){
  diag_R[i]<- sum(R[i,]== -1)
}

diag(R)<- diag_R
dim(R)

# SIMULATING NEW DATA:
{
  # taking true rho value as 0, independent case, generate spatial random effect psi
  n=50
  rho_true<- 0.75
  D_g_inv<-(1-rho_true)*diag(n)+rho_true*R
  sigma_true<- 1
  D<-solve(D_g_inv)
  set.seed(100)
  psi_true<- mvrnorm(1,rep(0,n),Sigma=D*sigma_true)
  x1<-rnorm(50)
  x<- cbind(1,x1)
  beta_true<- c(5,12)
  xbeta_true<- x%*% beta_true
  y<- xbeta_true+psi_true
  # generating response variable
  y<-(xbeta_true+psi_true)
  y<-round(y)
}

# PRE-MCMC:
{  
  # packages
  library(gmm)
  library(MASS)
  library(tidyverse)
  library(spdep)
  library(tidyverse)
  #install.packages("emplik")
  library(emplik)
  library(ape)
  
  # Functions used in the algorithm
  # Functions used in BEL leroux
  #estimating initial beta using el estimating equation taking psi=0,wi=(1/n), sigma_i^2=var(y)
  g1<-function(tet,x)
  {
    n<-dim(x)[1]
    beta1=(1/n)* sum(y-x%*%tet)
    beta2=(1/n)* (sum((y-x%*%tet)^2/var) -1)
    beta = cbind(beta1,beta2)
    return(beta)
  }
  
  dg<-function(tet,x)
  {
    n<-dim(x)[1]
    xx<-t(x)%*%x
    yx<-t(x)%*%y
    G<-matrix(c((1/n)*sum(-x),
                (2/n)*sum((t(x)%*%((x%*%tet)/var)-t(x)%*%(y/var)))),nrow=2,ncol=1) 
    return(G)
  }
  
  # calculating updated Beta using the old psi and wi_mu from last step to find wi*
  
  #var_orig<-as.numeric(var(y- x%*%beta-psi))
  #var <- var_orig
  
  
  #functions to calculate metropolis hastings ratio for accepting the proposed parameters
  # log(p(beta|Data,w,tau))
  # b_mean = prior mean
  target_beta<-function(b,w,b_mean,g,tau){
    sum(log(w))-0.5*(t(b-b_mean)%*%(b-b_mean))*g*tau
  }
  
  # log(p(psi|Data,w,tau))
  target_psi<-function(w,psi,D,tau){
    sum(log(w)) -0.5*t(psi)%*%D%*%psi*tau
  }
  
  # log(p(beta|tau)) *prior*
  beta_giv_tau<-function(b,b_mean,tau,g){
    btb<-as.vector(t(b-b_mean)%*%(b-b_mean))
    btb_tau<-btb*g*tau
    (-0.5*btb_tau)
  }
  
  psi_giv_tau<-function(psi,tau,D){-0.5*((psi%*%D)%*%psi)*tau}
  
  tau_prior<-function(tau,alpha_1,alpha_2){(1+alpha_1)*log(tau)-(alpha_2/tau)}
  
  target_tau<-function(tau,beta,beta_mean,alpha_1,alpha_2,psi,D){
    psi_giv_tau(psi,tau,D)+ beta_giv_tau(beta,beta_mean,tau,g)+tau_prior(tau,alpha_1,alpha_2)
  }
  
  # initialisation
  
  # initial tau
  n<- length(y) # no. of observations
  p<- dim(x)[2] # no. of covariates
  alpha_1<-1 # hyperparamter for tau prior
  alpha_2<-1 # hyperparamter for tau prior
  tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv or sigma^2, variance parameter
  tau_init<- 1/tau_inv_init
  
  g<- 10# G prior evaluated at 10 for regression coefficients' prior (Zellner prior)
  
  #initial Beta
  prior_mean_beta<- rep(0,p) # p is the number of regression parameters, in case of one covariate, p=2
  #beta_init<- rnorm(1,prior_mean_beta, (1/g)*tau_inv_init*diag(p))
  beta_init<- rnorm(2,prior_mean_beta, (1/g)*tau_inv_init)
  
  # MVN prior for Beta (check)
  #beta_init<-beta_true
  #beta_init<- model_fin$Beta[,1000]
  wi_init<- 1/length(y) # y be the response variable from the data
  # initial psi
  #psi_init<- psi_true
  psi_init <- rep(0,n) 
  var<- as.numeric(var(y- x%*%beta_init))
  ## MH algorithm for BEL leroux
  #niter= 100#  no. of iterations
  
  # calculating MELE of Beta, beta_mele
  wi=wi_init
  
  beta_mele<- unname(gel(g = g1, x = x, tet0= beta_init, gradv = dg)$coefficients) # caclulating MELE of Beta using gmm package
  
  # starting value of mu
  
  mu_init<- x%*% beta_mele + psi_init
  beta<-beta_mele
  
  # setting initial counter
  n.psi<- 0
  n.beta<- 0
  n.tau<- 0
  n.rho<-0
  # storing posterior samples
  
  wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
  wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 
  
  wi<-wi_mu
  
}

psi_init <- rnorm(50)
niter = 2000
rho_init = 0.2
sd_psi=3
sd_beta=3.5
sd_tau=0.99
sd_rho <- 5
#y,x,n,p,var,rho_init,niter,beta_init, psi_init, tau_init,R, wi, sd_psi, sd_beta, sd_tau

psi_sample<- matrix(nrow= n, ncol=niter)
beta_sample <- matrix(nrow= p, ncol=niter)
tau_sample <- c()
rho_sample<- c()
psi_sample[,1]<- psi_init
beta_sample[,1]<-beta
tau_sample[1]<-tau_init
rho_sample[1]<- rho_init
beta<-beta_sample[,1]
psi<-psi_sample[,1]
tau<-tau_sample[1]
rho<- rho_sample[1]
# setting a loop for iteration
i = 2
for(i in 2: niter){
  proposal.mean.psi<- psi_sample[,i-1]
  proposal.sd.psi<-sd_psi*diag(n)
  psi_star<-mvrnorm(1,proposal.mean.psi,proposal.sd.psi)
  
  # calculating updated Beta using the new psi and wi_mu from last step to find wi*
  #functions for estimating beta taking psi and wi from each step 
  g2_star<-function(tet,x, phi=psi_star)
  {
    beta1=sum(wi*(y-x%*%tet-phi))
    beta2=(sum(wi*((y-x%*%tet-phi)^2/var)) -1)
    beta = cbind(beta1,beta2)
    return(beta)
  }
  
  dg2_star<-function(tet,x,phi=psi_star)
  {
    n<-length(y)
    xx<-t(x)%*%x
    yx<-t(x)%*%y
    if(length(wi)==1){
      wi<-rep(wi,n)
    }
    G<-matrix(c(sum(-t(x)%*%wi)),
              (-2*sum((t(x)%*%(wi/var))%*%(t(wi)%*%(y-x%*%tet-phi)))),nrow=2,ncol=1)
    return(G)
  }
  beta<- beta_sample[,i-1]
  mu_new<-x%*%beta+ psi_star
  wi_star<-el.test(y-mu_new, 0)$wts
  wi_star<-wi_star/sum(wi_star)
  
  mu_old<-x%*%beta + psi
  wi_orig<-el.test(y-mu_old, 0)$wts
  wi_orig<-wi_orig/sum(wi_orig)
  
  # checking the constraint 14
  if(all(wi_star>0) & (sum(wi_star)-1) < 0.0001)
  {
    # perform a MH step for psi* in block k
    D_g_inv<-(1-rho)*diag(n)+rho*R # specify rho and neighbourhood matrix R beforehand
    pdr_psi<- target_psi(w=wi_star, psi=psi_star, D=D_g_inv, tau=tau)-
      target_psi(w=wi_orig, psi=psi, D=D_g_inv, tau=tau)#posterior density ratio for step 2, sampling psi
    #print(pdr_psi)
    if(rexp(1) > -pdr_psi){
      psi<-psi_star
      wi<- wi_star
      n.psi<-n.psi+1
    }
    
    else{
      psi<-psi
      wi<-wi
    }
  }
  
  psi_sample[,i]<-psi
  
  # # rho sampling
  # D_g_inv<-(1-rho)*diag(n)+rho*R
  # rho_proposed<-runif(1,0,1)
  # D_g_inv_proposed<-(1-rho_proposed)*diag(n)+rho_proposed*R
  # # sum(log(wi)) -0.5*t(psi)%*%D_g_inv%*%psi*tau
  # # sum(log(wi)) -0.5*t(psi)%*%D_g_inv_proposed%*%psi*tau
  # 
  # pdr_rho<- target_psi(w=wi, psi=psi, D=D_g_inv_proposed, tau=tau)-
  #   target_psi(w=wi, psi=psi, D=D_g_inv, tau=tau)
  # if(rexp(1) > -pdr_rho){
  #   rho<- rho_proposed
  #   n.rho<-n.rho+1
  # }else{
  #   rho<-rho
  # }
  
  {
    # rho sampling
    D_g_inv<-(1-rho)*diag(n)+rho*R
    lgt_rho_proposed<-rnorm(1, mean = logit(rho), sd=sd_rho)
    rho_proposed <- inv.logit(lgt_rho_proposed)
    D_g_inv_proposed<-(1-rho_proposed)*diag(n)+rho_proposed*R
    pdr_rho<- target_psi(w=wi, psi=psi, D=D_g_inv_proposed, tau=tau) - (log(rho) + log(1-rho)) -
      target_psi(w=wi, psi=psi, D=D_g_inv, tau=tau) + (log(rho_proposed) + log(1-rho_proposed))
    if(rexp(1) > -pdr_rho){
      rho<- rho_proposed
      n.rho<-n.rho+1
    }else{
      rho<-rho
    }
  }
  
  rho_sample[i]<- rho
  
  # Step 3 : sampling fixed effect beta
  
  # Step 3 : sampling fixed effect beta
  proposal.mean.beta<-unname(lm(y~x-1, weights = wi)$coefficients)
  proposal.var.beta<- sd_beta*diag(p) # specifying arbitrarily, large variance but not huge
  beta_proposed<- mvrnorm(1,proposal.mean.beta,proposal.var.beta)
  
  wi_beta<- el.test(y-x%*%beta_proposed-psi, 0)$wts
  wi_beta<- wi_beta/sum(wi_beta)
  wi_orig_2<-el.test(y-x%*%beta-psi,0)$wts
  wi_orig_2<-wi_orig_2/sum(wi_orig_2)
  
  
  if(all(wi_beta)>0 & (sum(wi_beta)-1)< 0.0001)
  {
    pdr_beta<-target_beta(beta_proposed,w=wi_beta,proposal.mean.beta,g=10,tau)-
      target_beta(beta,w=wi_orig_2,proposal.mean.beta,g=10,tau)# posterior density ratio for beta
    
    if(rexp(1) > -pdr_beta){
      wi<- wi_beta
      beta<-beta_proposed
      n.beta<-n.beta+1
    }
    else{ 
      wi<-wi
      beta<-beta
      n.beta<-n.beta}
  }
  beta_sample[,i]<-beta
  
  # step 4 : sampling precision parameter tau
  #mean_tau<-1/rgamma(1,alpha_1,alpha_2)
  #install.packages("truncnorm")
  #library(truncnorm)
  tau_proposed<- exp(rnorm(1,mean= log(tau),sd=sd_tau)) # random walk 
  
  
  pdr_tau<- target_tau(tau_proposed,beta_sample[,i],beta_mele,1,1,psi_sample[,i],D_g_inv)-
    target_tau(tau,beta_sample[,i-1],beta_mele,1,1,psi_sample[,i-1],D_g_inv) # posterior density ratio
  
  
  if(rexp(1)> -(pdr_tau+log(tau_proposed)-log(tau))) # if(rexp(1)> pdr_tau)
  {
    tau<-tau_proposed
    n.tau<-n.tau+1
  }  
  
  tau_sample[i]<- tau
}

par(mfrow = c(2,2))
plot(logit(rho_sample),type = "l")
plot(density(logit(rho_sample)))
plot((rho_sample),type = "l")
plot(density((rho_sample)))

cor(apply(psi_sample,1,mean), psi_true)


library(coda)
1-coda::rejectionRate(as.mcmc(t(psi_sample)))[1]
1-coda::rejectionRate(as.mcmc(logit(rho_sample)))
