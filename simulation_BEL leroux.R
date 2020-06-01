.libPaths("c:/software/Rpackages")
setwd("C:/R dir")
## generate some covariates
library(emplik)
x1<- runif(2147)
#x2<-rbinom(2147,1,0.5)
psi_true<-rnorm(2147)
beta_true<- c(7,14)
x<- cbind(1,x1)
xbeta_true<- x%*% beta_true
hist(xbeta_true)
summary(xbeta_true)
y<-exp(xbeta_true+psi_true)
y<-round(y)
summary(y)
hist(y)
hist(log(y))



# reading spatial data
##############################
library(tidyverse)
#install.packages("spdep")
library(spdep)
aus_nb<-read.gal("SA2_2011_AUST_qv7.gal", override.id = TRUE)
class(aus_nb)
W<-nb2mat(aus_nb,style="B")
nblist<-nb2listw(aus_nb)
#creating symmetric neighbourhood matrix for BYM in CARBAYES
#rownames(W)<-c()
#ind <- upper.tri(W)
#W[ind] <- t(W)[ind] 
#isSymmetric(W)
ni<-rowSums(W) # no. of neighbours for each area
R<-diag(ni)
for(i in 1:nrow(R))
{
  R[i,which(W[i,]==1)]<- -1
}

# reading cancer data

lung<-read_csv("lung.data.csv")
cov<-read_csv("SA2_v2.csv")

# parametric leroux

library(CARBayes)
param_fit<- S.CARleroux(formula= y~x-1,family ="poisson", fix.rho = T,rho = 0.1, W=W, burnin=100,n.sample=1000, thin= 1 )
    
##
Beta_true<-c(5,10)
# prior sd 5

summary(output$Beta[1,1000:2000])

