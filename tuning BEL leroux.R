

# finding acceptance
t1<- Sys.time()
model_1_.00025_3.25_0.99<-BEL_leroux_new(y,x,n,p,var,rho=0.1,niter=100,beta_init, psi_init, tau_init,R,wi,
                                     sd_psi=0.00025, sd_beta=3.5, sd_tau=0.99)
#save(model_red_sed_noblk3, file= "lerouxBEL_noblk_3000it.RData")
t2<-Sys.time()
t2-t1
model_1_.00022_3.25_0.99$acceptance_rate
model_1_.00025_3.5_0.99$acceptance_rate
plot(1:100,model_1_.00025_3.5_0.99$Beta[2,], type="l")
model_1_.00022_3.22_0.99$acceptance_rate
model_1_.00024_3.22_0.99$acceptance_rate
model_1_.00025_3.22_0.99$acceptance_rate
model_1_.00025_3.22_0.99$acceptance_rate

