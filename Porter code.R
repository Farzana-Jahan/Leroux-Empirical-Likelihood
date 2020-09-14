## Functions to be used

#estimating initial beta using el estimating equation taking Yq=0
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
              (2/(n*var))*sum((t(x)%*%((x%*%tet))-yx))),nrow=2,ncol=1) 
  return(G)
}

#estimating equations for beta using Miyq

g2<-function(tet,x)
{
  wi<-el.test(y-x%*%tet-miyq,0)$wts
  wi<-wi/sum(wi)
  beta1=sum(wi*(y-x%*%tet-miyq))
  beta2=(sum(wi*((y-x%*%tet-miyq)^2/var)) -1)
  beta = cbind(beta1,beta2)
  return(beta)
}

dg2<-function(tet,x)
{
  n<-length(y)
  wi<-el.test(y-x%*%tet-miyq,0)$wts
  wi<-wi/sum(wi)
  xx<-t(x)%*%x
  yx<-t(x)%*%y
  G<-matrix(c(sum(-t(x)%*%wi)),
            (-2*sum((t(x)%*%(wi/var))%*%(t(wi)%*%(y-x%*%tet-miyq)))),nrow=2,ncol=1)
  return(G)
}


#functions for metropolis hastings

target_yq<-function(a,w,MBM,tau){sum(log(w))-0.5*t(a)%*%MBM%*%a*tau}


target_beta<-function(b,w,b_mean,g,tau){
  sum(log(w))-0.5*(t(b-b_mean)%*%(b-b_mean))*g*tau
}

yq_giv_tau<-function(a,tau,MBM){as.numeric(-0.5*(a%*%MBM%*%a)*tau)}

tau_prior<-function(tau,alpha_1,alpha_2){(1+alpha_1)*log(tau)-(alpha_2/tau)}

beta_giv_tau<-function(b,b_mean,tau,g){
  btb<-as.vector(t(b-b_mean)%*%(b-b_mean))
  btb_tau<-btb*g*tau
  -0.5*btb_tau
}


target_tau<-function(tau,beta,beta_mean,alpha_1,alpha_2,psi,MBM,p,q){
 ((p+q)/2)*log(tau)+ yq_giv_tau(psi,tau,MBM)+ beta_giv_tau(beta,beta_mean,tau,g)+tau_prior(tau,alpha_1,alpha_2)
}

#finding mode of fixed effects
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#creating Moran's basis and other required measures
M_create=function(y,x,B, Bplus)
{
  n<-length(y)
  xx<-t(x)%*%x
  p=diag(1,n)-x%*%solve(xx)%*%t(x)
  pbp=t(p)%*%B%*%p
  q<-sum(eigen(pbp)$values>0)
  Re(eigen(pbp)$vectors[,1:q])
}


MBM_create=function(M,B, Bplus)
{
  t(M)%*%(B-Bplus)%*%M
}

BSHEL<-function(y,x,n,p,var,niter,beta_init, psi_init, tau_init,M,MBM, wi, sd_psi, sd_beta, sd_tau)
{
  # setting initial counter
  n.psi<- 0
  n.beta<- 0
  n.tau<- 0
  q=dim(MBM)[2]
  psi_sample<- matrix(0,nrow= q, ncol=niter)
  beta_sample<- matrix(0,nrow= p, ncol=niter)
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
    beta<- beta_sample[,i-1]
    
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
     # specify rho and neighbourhood matrix R beforehand
      pdr_psi<- as.numeric(target_yq(w=wi_star, a=psi_star, MBM, tau=tau)-
        target_yq(w=wi_orig, a=psi, MBM, tau=tau))#posterior density ratio for step 2, sampling psi
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
    
    wi_beta<- el.test(y-x%*%beta_proposed-M%*%psi, 0)$wts
    wi_beta<- wi_beta/sum(wi_beta)
    wi_orig_2<-el.test(y-x%*%beta-M%*%psi,0)$wts
    wi_orig_2<-wi_orig_2/sum(wi_orig_2)
    
    
    if(all(wi_beta)>0 & (sum(wi_beta)-1)< 0.0001)
    {
      pdr_beta<-target_beta(beta_proposed,w=wi_beta,proposal.mean.beta,g=10,tau)-
        target_beta(beta,w=wi_orig_2,proposal.mean.beta,g=10,tau)# posterior density ratio for beta
      if(rexp(1) > -pdr_beta){
        wi<- wi_beta
        beta<-beta_proposed
        n.beta<-n.beta+1
        var<- as.numeric(var(y- x%*%beta))
      }
      else{ 
        wi<-wi
        beta<-beta
        n.beta<-n.beta
        var<- var}
    }
    beta_sample[,i]<-beta
    
    # step 4 : sampling precision parameter tau
    #mean_tau<-1/rgamma(1,alpha_1,alpha_2)
    #install.packages("truncnorm")
    #library(truncnorm)
    tau_proposed<- exp(rnorm(1,mean= log(tau),sd=sd_tau)) # random walk 
    
    
    pdr_tau<- target_tau(tau_proposed,beta_sample[,i],proposal.mean.beta,1,1,psi_sample[,i],MBM,p,q)-
      target_tau(tau,beta_sample[,i-1],proposal.mean.beta,1,1,psi_sample[,i-1],MBM,p,q) # posterior density ratio
    
    
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
  output<- list(Beta= beta_sample, psi= psi_sample, tau= tau_sample, acceptance_rate=acceptance, 
                effective_sample_size=ess)
  output
}


# Data reading from North Carolina SIDS 

library(gmm)
library(MASS)
library(tidyverse)
library(spdep)
library(tidyverse)
#install.packages("emplik")
library(emplik)
library(ape)


#install.packages("emplik")
library(emplik)
# calculating MELE's of fixed effects beta using MELE by GMM (gel methods) 

library(gmm)
library(rgdal)
## from data
# from data set:
map <- readOGR(system.file("shapes/sids.shp", package = "spData")[1])
# View data contained within shapefile
data <- as.data.frame(map)
head(data)
data <- data[,-c(1:7, 15:22)]
head(data)
dat.1 <- data.frame(data[,c(1:4)], year = "1974-78")
names(dat.1) <- c("area.ID", "pop", "observed", "x", "year")
dat.2 <- data.frame(data[,c(1, 5:7)], year = "1979-84")
names(dat.2) <- names(dat.1)
data <- rbind(dat.1, dat.2)
rm(dat.1, dat.2)
head(data)
data$expected <- c(
  data$pop[1:100] * sum(data$observed[1:100]) / sum(data$pop[1:100]),
  data$pop[101:200] * sum(data$observed[101:200]) / sum(data$pop[101:200])
)
data$raw <- data$observed / data$expected
data_1<-data[data$year=="1974-78",]
head(data_1)
summary(data_1)
data_1$raw[data_1$raw==0]<- 0.01
data_1$x<- data_1$x/100

# creating neighbourhood matrix

sids_nb<-read.gal("C:/R dir/Leroux-Empirical-Likelihood-master/Leroux-Empirical-Likelihood-master/shapefile_SIDS.gal", override.id = TRUE)
#class(sids_nb)
W<-nb2mat(sids_nb,style="B")
nblist<-nb2listw(sids_nb)
#creating symmetric neighbourhood matrix for BYM in CARBAYES
rownames(W)<-c()
ind <- upper.tri(W)
W[ind] <- t(W)[ind] 
#isSymmetric(W)
ni<-rowSums(W)


# For BSHEL

y<- log(data_1$raw)
x<- cbind(1,data_1$x)
B<-W
B_plus<-diag(rowSums(B))

n<-length(y)
M=M_create(y,x,B,B_plus)
MBM=MBM_create(M,B,B_plus)

q=dim(MBM)[2]
psi_init <- rep(0,q) 

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
var<- as.numeric(var(y- x%*%beta_init))
# MVN prior for Beta (check)
#beta_init<-beta_true
#beta_init<- model_fin$Beta[,1000]
wi_init<- 1/length(y) # y be the response variable from the data
psi_init<- rep(0,q)



# calculating MELE of Beta, beta_mele
wi=wi_init

beta_mele<- unname(gel(g = g1, x = x, tet0= beta_init, gradv = dg)$coefficients) # caclulating MELE of Beta using gmm package

# starting value of mu

mu_init<- x%*% beta_mele + M%*%psi_init
beta_init<-beta_mele
var<- as.numeric(var(y- x%*%beta_init))

# storing posterior samples

wi_mu<- el.test(y-mu_init,0)$wts # computing el weights using emplik package
wi_mu<-wi_mu/sum(wi_mu) # sum(wi) = 1 and wi>0 constraints 

wi<-wi_mu

## Running porter model (BSHEL) for this data

porter_1k<- BSHEL(y,x,n,p,var,niter=10000,beta_init, psi_init, tau_init,M,MBM, wi, sd_psi=0.01, sd_beta=0.0000025, sd_tau=1.5)
