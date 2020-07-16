# creating a autoregression matrix randomly to model BYM and leroux 

# creating autoregression matrix randomly for a set of 50 small areas

R<- matrix(sample(-1:0,2500,replace=T),nrow=50, ncol=50)
diag_R<-c()
for(i in 1:50){
  diag_R[i]<- sum(R[i,]== -1)
}

diag(R)<- diag_R
dim(R)

# taking true rho value as 0, independent case, generate spatial random effect psi
n=50
rho_true<- 0
D_g_inv<-(1-rho_true)*diag(n)+rho_true*R
sigma_true<- 0.2
D<-solve(D_g_inv)
set.seed(100)
psi_true<- mvrnorm(1,rep(0,n),Sigma=D*sigma_true)
x1<-rnorm(50)
x<- cbind(1,x1)
beta_true<- c(5,12)
xbeta_true<- x%*% beta_true
y<- xbeta_true+psi_true
# generating response variable
y<-exp(xbeta_true+psi_true)
y<-round(y)

# taking true rho value as 0.75, moderate to high autocorrelation case, generate spatial random effect psi
n=50
rho_true<- 0.75
D_g_inv<-(1-rho_true)*diag(n)+rho_true*R
sigma_true<-1
D<-solve(D_g_inv)
# install.packages("corpcor")
library(corpcor)
set.seed(100)
psi_true<- mvrnorm(1,rep(0,n),Sigma= make.positive.definite(D*sigma_true))
x1<-rnorm(50)
x<- cbind(1,x1)
beta_true<- c(5,12)
xbeta_true<- x%*% beta_true
y<- xbeta_true+psi_true
y<-round(y)



