library(MASS)

R<- matrix(sample(-1:0,2500,replace=T),nrow=50, ncol=50)
diag_R<-c()
for(i in 1:50){
  diag_R[i]<- sum(R[i,]== -1)
}

diag(R)<- diag_R
dim(R)
R[1:10,1:10]
ind <- upper.tri(R)
R[ind] <- t(R)[ind] 
isSymmetric(R) # have to make R symmetric

# creating B and B plus for POrter, B in porter is equivalent to R is leroux
B_plus<-diag(diag_R)
B_plus[1:10,1:10]
B<-R
B[B==-1]<- 1
B[1:10,1:10]

# data generation
# taking true rho value as 0, independent case, generate spatial random effect psi
library(corpcor)
n=50
rho_true<- 0.7
D_g_inv<-(1-rho_true)*diag(n)+rho_true*R
sigma_true<- 0.2
tau_true<- 1/sigma_true
D<-ginv(tau_true*D_g_inv)
set.seed(100)
psi_true<- mvrnorm(1,rep(0,n),Sigma=make.positive.definite(D*sigma_true))
x1<-rnorm(50)
x<- cbind(1,x1)
xx<- t(x)%*%x
beta_true<- c(5,12)
xbeta_true<- x%*% beta_true
y<- xbeta_true+psi_true
# generating response variable
#y<-exp(xbeta_true+psi_true)
#y<-exp(xbeta_true+psi_true)
y<-round(y)

# creating M and MBM
# functions to create M and MBM
#creating Moran's basis and other required measures
M_create=function(y,x,B, Bplus)
{
  n<-length(y)
  p=diag(1,n)-x%*%solve(xx)%*%t(x)
  pbp=t(p)%*%B%*%p
  q<-sum(eigen(pbp)$values>0)
  Re(eigen(pbp)$vectors[,1:q])
}


MBM_create=function(M,B, Bplus)
{
  t(M)%*%(B-Bplus)%*%M
}
M=M_create(y,x,B,B_plus)
MBM=MBM_create(M,B,B_plus)

#step 2 to 4 for Porter paper
library(gmm)
library(MASS)
library(tidyverse)
library(spdep)
library(tidyverse)
#install.packages("emplik")
library(emplik)
library(ape)

# Functions used in the algorithm

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

#functions to calculate metropolis hastings ratio for accepting the proposed parameters
# log(p(beta|Data,w,tau))
# b_mean = mean of the proposal dist for beta, taken as the WLS estimates of beta from the data
target_beta<-function(b,w,b_mean,g,tau){
  sum(log(w))-0.5*(t(b-b_mean)%*%(b-b_mean))*g*tau
}

# log(p(psi|Data,w,tau))
target_psi<-function(w,psi,MBM,tau){
  sum(log(w)) -0.5*t(psi)%*%MBM%*%psi*tau
}

# log(p(beta|tau)) *prior*
beta_giv_tau<-function(b,b_mean,tau,g){
  btb<-as.vector(t(b-b_mean)%*%(b-b_mean))
  btb_tau<-btb*g*tau
  (-0.5*btb_tau)
}

psi_giv_tau<-function(psi,tau,MBM){-0.5*((psi%*%MBM)%*%psi)*tau}

tau_prior<-function(tau,alpha_1,alpha_2){(1+alpha_1)*log(tau)-(alpha_2/tau)}

target_tau<-function(tau,beta,beta_mean,alpha_1,alpha_2,psi,MBM){
  psi_giv_tau(psi,tau,MBM)+ beta_giv_tau(beta,beta_mean,tau,g)+tau_prior(tau,alpha_1,alpha_2)
}

# functions to create M and MBM
#creating Moran's basis and other required measures
M_create=function(y,x,B, Bplus)
{
  n<-length(y)
  p=diag(1,n)-x%*%solve(xx)%*%t(x)
  pbp=t(p)%*%B%*%p
  q<-sum(eigen(pbp)$values>0)
  Re(eigen(pbp)$vectors[,1:q])
}


MBM_create=function(M,B, Bplus)
{
  t(M)%*%(B-Bplus)%*%M
}
M=M_create(y,x,B,B_plus)
MBM=MBM_create(M,B,B_plus)


library(emplik)
# calculating MELE's of fixed effects beta using MELE by GMM (gel methods) 

library(gmm)

# initial tau
n<- length(y) # no. of observations
p<- dim(x)[2] # no. of covariates
alpha_1<-1 # hyperparamter for tau prior
alpha_2<-1 # hyperparamter for tau prior
tau_inv_init<- rgamma(1,alpha_1,alpha_2) # using IG prior(1,1) for tau_inv
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
q=dim(MBM)[2]
psi_init <- rep(0,q) ##BUG##: remove line when beta-proposal issue fixed.
var<- as.numeric(var(y- x%*%beta_init))
## MH algorithm for BEL leroux
#niter= 100#  no. of iterations

# calculating MELE of Beta, beta_mele
wi=wi_init

beta_mele<- unname(gel(g = g1, x = x, tet0= beta_init, gradv = dg)$coefficients) # caclulating MELE of Beta using gmm package

# starting value of mu

mu_init<- x%*% beta_mele + M%*%psi_init
beta_init<-beta_mele


# storing posterior samples

wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

wi<-wi_mu


BSHEL<-function(y,x,n,p,var,niter,beta_init, psi_init, tau_init,M,MBM, wi, sd_psi, sd_beta, sd_tau)
{
  # setting initial counter
  n.psi<- 0
  n.beta<- 0
  n.tau<- 0
  q=dim(MBM)[2]
  psi_sample<- matrix(nrow= q, ncol=niter)
  beta_sample<- matrix(nrow= p, ncol=niter)
  tau_sample<- c()
  psi_sample[,1]<- psi_init
  beta_sample[,1]<-beta_init
  tau_sample[1]<-tau_init
  beta<-beta_sample[,1]
  psi<-psi_sample[,1]
  tau<-tau_sample[1]
  # setting a loop for iteration
  for(i in 2: niter){
    
    proposal.mean.psi<- psi_sample[,i-1]
    proposal.sd.psi<-sd_psi*diag(q)
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
    mu_new<-x%*%beta+ M%*%psi_star
    wi_star<-el.test(y-mu_new, 0)$wts
    wi_star<-wi_star/sum(wi_star)
    
    mu_old<-x%*%beta + M%*%psi
    wi_orig<-el.test(y-mu_old, 0)$wts
    wi_orig<-wi_orig/sum(wi_orig)
    
    # checking the constraint 14
    if(all(wi_star>0) & (sum(wi_star)-1) < 0.0001)
    {
      # perform a MH step for psi* in block k
      pdr_psi<- target_psi(w=wi_star, psi=psi_star, MBM, tau=tau)-
        target_psi(w=wi_orig, psi=psi, MBM, tau=tau)#posterior density ratio for step 2, sampling psi
      if(rexp(1) > - pdr_psi){
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
    
    
    pdr_tau<- target_tau(tau_proposed,beta_sample[,i],beta_mele,1,1,psi_sample[,i],MBM)-
      target_tau(tau,beta_sample[,i-1],beta_mele,1,1,psi_sample[,i-1],MBM) # posterior density ratio
    
    
    if(rexp(1)> -(pdr_tau+log(tau_proposed)-log(tau))) # if(rexp(1)> pdr_tau)
    {
      tau<-tau_proposed
      n.tau<-n.tau+1
    }  
    
    tau_sample[i]<- tau
  }
  
  acceptance<- data.frame(parameter= c("beta", "psi", "tau"), 
                          acceptance_rate= c((n.beta/niter)*100, (n.psi/niter)*100,
                                             (n.tau/niter)*100))
  n.psi<-n.psi
  # effective sample size 
  # autocorrelation
  acf_psi<-c()
  for(j in 1:q)
  {
    acf<-acf(psi_sample[j,],plot=F)[[1]]
    n<-length(acf)
    acf_psi[j]<-sum(acf[2:q])
  }
  acf_beta0<-acf(beta_sample[1,],plot=F)[[1]]
  acf_beta1<-acf(beta_sample[2,],plot=F)[[1]]
  acf_tau<-acf(tau_sample,plot=F)[[1]]
  # ess
  ess_psi<-c()
  for(j in 1:q)
  {
    ess_psi[j]<- niter/ (1+2*acf_psi[j])
  }
  ess_beta0<- n.beta/ (1+2*sum(acf_beta0))
  ess_beta1<- n.beta/ (1+2*sum(acf_beta1))
  ess_tau<-n.tau/(1+2*sum(acf_tau))
  ess<- list(ess_psi, ess_beta0, ess_beta1, ess_tau)
  # tau, beta, best and worst psi
  output<- list(Beta= beta_sample, psi= psi_sample, tau= tau_sample, acceptance_rate=acceptance, 
                effective_sample_size=ess)
  output
}

# run the model:

model_porter_1000<-BSHEL(y,x,n,p,var,niter=1000,beta_init, psi_init, tau_init,M,MBM, wi, sd_psi=6, sd_beta=0.000005, sd_tau=1.6)
model_porter_1000$acceptance_rate



 
