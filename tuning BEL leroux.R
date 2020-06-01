

# finding reasonable acceptance for beta, psi and tau in BEL leroux
# sd psi = 0.0025, sd_beta= 3.5, sd_tau= 0.99
t1<- Sys.time()
model_1_.00025_3.5_0.99<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=500,beta_init, psi_init, tau_init,R,wi,
                                     sd_psi=0.00025, sd_beta=3.5, sd_tau=0.99)
t2<-Sys.time()
t2-t1
model_1_.00025_3.5_0.99$acceptance_rate
#  trace plots for beta
plot(1:500,model_1_.00025_3.5_0.99$Beta[1,], type="l")
plot(1:500,model_1_.00025_3.5_0.99$Beta[2,], type="l")

# sd psi = 0.0025, sd_beta= 3.25, sd_tau= 0.85
t1<- Sys.time()
model_1_.00025_3.25_0.85<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=500,beta_init, psi_init, tau_init,R,wi,
                                     sd_psi=0.00025, sd_beta=3.25, sd_tau=0.85)
t2<-Sys.time()
t2-t1
model_1_.00025_3.25_0.85$acceptance_rate
#  trace plots for beta
plot(1:500,model_1_.00025_3.25_0.85$Beta[1,], type="l")
plot(1:500,model_1_.00025_3.25_0.85$Beta[2,], type="l")

# sd psi = 0.0024, sd_beta= 3, sd_tau= 0.75
t1<- Sys.time()
model_1_.00025_3_0.75<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=500,beta_init, psi_init, tau_init,R,wi,
                                     sd_psi=0.00025, sd_beta=3, sd_tau=0.75)
t2<-Sys.time()
t2-t1
model_1_.00025_3_0.75$acceptance_rate
#  trace plots for beta
plot(1:500,model_1_.00025_3_0.75$Beta[1,], type="l")
plot(1:500,model_1_.00025_3_0.75$Beta[2,], type="l")
# sd psi = 0.0024, sd_beta= 3, sd_tau= 0.65
t1<- Sys.time()
model_1_.00024_2_0.65<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=500,beta_init, psi_init, tau_init,R,wi,
                                     sd_psi=0.00024, sd_beta=2, sd_tau=0.65)
t2<-Sys.time()
t2-t1
model_1_.00024_2_0.65$acceptance_rate
#  trace plots for beta
plot(1:500,model_1_.00024_2_0.65$Beta[1,], type="l")
plot(1:500,model_1_.00024_2_0.65$Beta[2,], type="l")


