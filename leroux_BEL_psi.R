
# packages
library(gmm)
library(MASS)
library(tidyverse)
library(spdep)
library(tidyverse)
#install.packages("emplik")
library(emplik)

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
tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
tau_init<- 1/tau_inv_init
#tau_init<-model_fin_new2$tau[5000]
#tau_init <- 1/var ##BUG##: remove line when beta-proposal issue fixed.
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
psi_init<- psi_true
#psi_init <- 0 ##BUG##: remove line when beta-proposal issue fixed.
var<-as.numeric(var(psi_true))
#var<- as.numeric(var(y- x%*%beta_init))
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
# storing posterior samples

wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

wi<-wi_mu

BEL_leroux_new_psi<-function(y,x,n,p,var,rho,niter,beta_init, psi_init, tau_init,R, wi, sd_psi, sd_beta, sd_tau)
{
  
  psi_sample<- matrix(nrow= n, ncol=niter)
  beta_sample <- matrix(nrow= p, ncol=niter)
  tau_sample <- c()
  psi_sample[,1]<- psi_init
  beta_sample[,1]<-beta
  tau_sample[1]<-tau_init
  beta<-beta_sample[,1]
  psi<-psi_sample[,1]
  tau<-tau_sample[1]
  # Creating matrices for proposal kernels:
  proposal.var.beta<- sd_beta*diag(p) # specifying arbitrarily, large variance but not huge
  #proposal.sd.psi<-sd_psi*diag(n)
  t1 <- proc.time()
  # setting a loop for iteration
  for(i in 2: niter){
    
    proposal.mean.psi<- rep(0,n)
    beta<- beta_sample[,i-1]
    psi_star<-rnorm(length(proposal.mean.psi),proposal.mean.psi,sd_psi)
 
    mu_new<-x%*%beta + psi_star
    wi_star<-el.test(y-mu_new, 0)$wts
    wi_star<-wi_star/sum(wi_star)
    
    
    # calculating updated Beta using the old psi and wi_mu from last step to find wi*
 
    mu_old<-x%*%beta + psi
    wi_orig<-el.test(y-mu_old, 0)$wts
    wi_orig<-wi_orig/sum(wi_orig)
    
    # # checking the constraint 14
    if(all(wi_star>0) & (sum(wi_star)-1) < 0.0001)
    {
    # perform a MH step for psi* 
    D_g_inv<-(1-rho)*diag(n)+rho*R # specify rho and neighbourhood matrix R beforehand
    pdr_psi<- target_psi(w=wi_star, psi=psi_star, D=D_g_inv, tau=tau)/
    target_psi(w=wi_orig, psi=psi, D=D_g_inv, tau=tau)#posterior density ratio for step 2, sampling psi
    if(rexp(1) > -pdr_psi){
    psi<-psi_star
    wi<- wi_star
    n.psi<-n.psi+1
    }
       
    else{
    psi<-psi
    wi<-wi
    n.psi<-n.psi
    }
    }
   psi_sample[,i]<-psi
    
    # Step 3 : sampling fixed effect beta
    #var= as.numeric(var(y- x%*%beta-psi_sample[,i]))
   g2<-function(tet,x, phi=psi)
   {
     beta1=sum(wi*(y-x%*%tet-phi))
     beta2=(sum(wi*((y-x%*%tet-phi)^2/var)) -1)
     beta = cbind(beta1,beta2)
     return(beta)
   }
   
   dg2<-function(tet,x,phi=psi)
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
   
   
    #proposal.mean.beta<- unname(gel(g2,x,beta,gradv = dg2)$coefficients)
   proposal.mean.beta<-unname(lm(y~x-1, weights = wi)$coefficients)
   # beta_proposed<- mvrnorm(1,beta,proposal.var.beta)
    beta_proposed<- rnorm(2,proposal.mean.beta,sd_beta)
   # print(beta_proposed)
    wi_beta<- el.test(y-x%*%beta_proposed-psi,0)$wts
    wi_beta<- wi_beta/sum(wi_beta)
    wi_orig_2<-el.test(y-x%*%beta-psi,0)$wts
    wi_orig_2<-wi_orig_2/sum(wi_orig_2)
    
    
    if(all(wi_beta)>0 & (sum(wi_beta)-1)< 0.0001)
    {
      pdr_beta<-target_beta(beta_proposed,w=wi_beta,proposal.mean.beta,g=10,tau)-
        target_beta(beta,w=wi_orig_2,proposal.mean.beta,g=10,tau)# posterior density ratio for beta
      #print(pdr_beta)
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
    mean_tau<-1/rgamma(1,alpha_1,alpha_2)
    #install.packages("truncnorm")
    library(truncnorm)
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
  t2 <- proc.time()
  t2[3] - t1[3]
  (t2[3] - t1[3])/(niter-1)
  acceptance<- data.frame(parameter= c("beta", "psi", "tau"), 
                          acceptance_rate= c((n.beta/niter)*100, (n.psi/niter)*100,
                                             (n.tau/niter)*100))
  # effective sample size 
  # autocorrelation
   acf_psi<-c()
  for(j in 1:2147)
  {
   acf<-acf(psi_sample[j,],plot=F)[[1]]
  n<-length(acf)
   acf_psi[j]<-sum(acf[2:n])
  }
  acf_beta0<-acf(beta_sample[1,],plot=F)[[1]]
  acf_beta1<-acf(beta_sample[2,],plot=F)[[1]]
  acf_tau<-acf(tau_sample,plot=F)[[1]]
  # ess
  ess_psi<-c()
  for(j in 1:2147)
  {
   ess_psi[j]<- niter/ (1+2*acf_psi[j])
  }
  ess_beta0<- n.psi/ (1+2*sum(acf_beta0))
  ess_beta1<- n.beta/ (1+2*sum(acf_beta1))
  ess_tau<-n.tau/(1+2*sum(acf_tau))
  ess<- list(ess_psi,ess_beta0, ess_beta1, ess_tau)
  #ess<-  list(ess_beta0, ess_beta1, ess_tau)
  # tau, beta, best and worst psi
  output<- list(Beta= beta_sample,psi= psi_sample, 
                tau= tau_sample, acceptance_rate=acceptance, 
                effective_sample_size=ess)
  output
}
