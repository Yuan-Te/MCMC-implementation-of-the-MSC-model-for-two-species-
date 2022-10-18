#Assignment 3
HC.SitesDiffs <- read.delim(file.choose())

L = 1000
dim(HC.SitesDiffs)
HC_1000sites <- HC.SitesDiffs[1:L,]
x <- c(HC_1000sites[,2])
n <- c(HC_1000sites[,1])



# tau = species divergence time
# theta = population size parameter
# N = number of samples including initial state
# x = number of observed differences between two sequences
# n = sequence length
# w = window size
# mu = prior mean of distance

lnPf = function(tau, theta, time, mu_tau, mu_theta, x, n){
  sigma = c()
  d=2*(tau+theta) #tau will renew
  p = 3/4 - (3/4)*exp(-(4*d)/3)
  sigma = (log(2/theta) - 2*time/theta + x*log(p) + (n-x)*log(1-p))
  lnP = -tau/mu_tau - theta/mu_theta + sum(sigma)
  return(lnP) #lnP is log(posterior)
}

lnposterior <- lnPf(tau = 0.01, theta = 0.001, time = 0.001, mu_tau = 0.005, mu_theta = 0.001, x, n)
lnposterior



msc = function(tau, theta, time, N, x, n, w_tau, w_theta, w_time, mu_tau, mu_theta){
  sample_tau = numeric(N+1)
  sample_theta = numeric(N+1)
  sample_t = matrix(nrow = N, ncol = L)
  sample_tau[1] = tau
  sample_theta[1] = theta
  sample_t[1,] = time
  acceptance_tau = 0
  acceptance_theta = 0
  acceptance_time = 0
  
  lnp = lnPf(tau, theta, time, mu_tau, mu_theta, x, n);

  for (i in 1:N) {
        
    # Step 2a (changing tau)
    tau_new = tau - w_tau/2 + w_tau*runif(1);
    if (tau_new < 0) tau_new = -tau_new;
    
    lnp_new = lnPf(tau_new, theta, time, mu_tau, mu_theta, x, n);
    logratio = lnp_new - lnp;
    if(logratio>=0 || runif(1)<exp(logratio)) {
        tau = tau_new;
        lnp = lnp_new;
        acceptance_tau = acceptance_tau+1
    }
    else {
      ;
    }
    sample_tau[i+1] = tau;

    # Step 2b (changing theta)
    theta_new = theta - w_theta/2 + w_theta*runif(1);
    if (theta_new < 0) theta_new = -theta_new;
    
    lnp_new = lnPf(tau, theta_new, time, mu_tau, mu_theta, x, n);
    logratio = lnp_new - lnp;
    
    if(logratio>=0 || runif(1)<exp(logratio)) {
        theta = theta_new;
        lnp = lnp_new;
        acceptance_theta = acceptance_theta+1
    }
    else {
      ;
    }
    sample_theta[i+1] = theta;

    # Step 2c (a loop over locus j changing the 1000 coalescent times)
    for(j in 1:L){
        t_new = time - w_time/2 + w_time*runif(1);
        if (t_new < 0) t_new = -t_new;
        lnp_new = lnPf(tau, theta, t_new, mu_tau, mu_theta, x, n);
        logratio = lnp_new - lnp;
        if(logratio>=0 || runif(1)<exp(logratio)) {
            time = t_new;
            lnp = lnp_new;
            acceptance_time = acceptance_time+1
        }
        sample_t[i,j] = time
    }


  }
  acceptance_tau = acceptance_tau/(N - 1)
  acceptance_theta = acceptance_theta/(N - 1)
  acceptance_time = acceptance_time/(L*(N-1))
  return(list(sample_tau, sample_theta, sample_t, acceptance_tau, acceptance_theta, acceptance_time))
}


#Window size
run5 = msc(tau=0.01,  theta=0.001,  time=0.001, N=1000, x, n, w_tau=0.0007, w_theta=0.00002, w_time=0.000001, mu_tau=0.005, mu_theta=0.001)
c(run5[[4]],run5[[5]],run5[[6]])
plot(run5[[1]], type="l", col="black", xlab="Iterations", ylab="Tau and Theta", ylim=c(0, 0.01), xlim=c(0, 1000))
lines(run5[[2]], col="red")

sample_tau = run5[[1]][100:length(run5[[1]])];
summary(sample_tau)  #posterior mean for 1st MCMC, dropping 100 values as burn-in
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_theta = run5[[2]][400:length(run5[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 400 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

#Efficiency
sample_tau.acf = acf(sample_tau)   # The correlation is lost after a lag of k=10
sample_tau.acf$acf 
eff = function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))
eff(sample_tau.acf)

sample_theta.acf = acf(sample_theta)   # The correlation is lost after a lag of k=10
sample_theta.acf$acf 
eff = function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))
eff(sample_theta.acf)

#Easy way for CI
sample_tau <- data.frame(sample_tau)
l.model <- lm(sample_tau~1,sample_tau)
confint(l.model, level=0.95)

sample_theta <- data.frame(sample_theta)
l.model <- lm(sample_theta~1,sample_theta)
confint(l.model, level=0.95)




#Modify tau, theta, time
run1 = msc(tau=0.01, theta=0.001, time=0.001, N=1000, x, n, w_tau=0.01, w_theta=0.002, w_time=0.002, mu_tau=0.005, mu_theta=0.001)
run2 = msc(tau=0.1,  theta=0.01,  time=0.01, N=1000, x, n, w_tau=0.01, w_theta=0.002, w_time=0.002, mu_tau=0.005, mu_theta=0.001)
run3 = msc(tau=0.3,  theta=0.05,  time=0.1, N=1000, x, n, w_tau=0.01, w_theta=0.002, w_time=0.002, mu_tau=0.005, mu_theta=0.001)
run4 = msc(tau=0.5,  theta=0.1,  time=4, N=1000, x, n, w_tau=0.01, w_theta=0.002, w_time=0.0002, mu_tau=0.005, mu_theta=0.001)
c(run1[[4]], run2[[4]], run3[[4]], run4[[4]]) #acceptance_tau
c(run1[[5]], run2[[5]], run3[[5]], run4[[5]]) #acceptance_theta
c(run1[[6]], run2[[6]], run3[[6]], run4[[6]]) #acceptance_time

#Tau
plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="Tau", ylim=c(0, 0.5), xlim=c(0, 1000))
lines(run2[[1]], col="red")
lines(run3[[1]], col="blue")
lines(run4[[1]], col="green")
#burn-in: 400

#Theta
plot(run1[[2]], type="l", col="black", xlab="Iterations", ylab="Theta", ylim=c(0, 0.2), xlim=c(0, 1000))
lines(run2[[2]], col="red")
lines(run3[[2]], col="blue")
lines(run4[[2]], col="green")
#burn-in: 800

#Impact of the priors on the posterior_tau
sample_tau = run1[[1]][400:length(run4[[1]])];
summary(sample_tau)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run2[[1]][400:length(run4[[1]])];
summary(sample_tau)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run3[[1]][400:length(run4[[1]])];
summary(sample_tau)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run4[[1]][400:length(run4[[1]])];
summary(sample_tau)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

#Impact of the priors on the posterior_theta
sample_theta = run1[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run2[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run3[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run4[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta



#Modify mu
run1 = msc(tau=0.01, theta=0.001, time=0.001, N=1000, x, n, w_tau=0.0007, w_theta=0.00002, w_time=0.000001, mu_tau=0.005, mu_theta=0.001)
run2 = msc(tau=0.01, theta=0.001, time=0.001, N=1000, x, n, w_tau=0.0007, w_theta=0.00002, w_time=0.000001, mu_tau=0.02, mu_theta=0.004)
run3 = msc(tau=0.01, theta=0.001, time=0.001, N=1000, x, n, w_tau=0.0007, w_theta=0.00002, w_time=0.000001, mu_tau=0.5, mu_theta=0.1)
run4 = msc(tau=0.01, theta=0.001, time=0.001, N=1000, x, n, w_tau=0.0007, w_theta=0.00002, w_time=0.000001, mu_tau=5, mu_theta=1)

c(run1[[4]], run2[[4]], run3[[4]], run4[[4]]) #acceptance_tau
c(run1[[5]], run2[[5]], run3[[5]], run4[[5]]) #acceptance_theta
c(run1[[6]], run2[[6]], run3[[6]], run4[[6]]) #acceptance_time

#Tau
plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="Tau", ylim=c(0, 0.01), xlim=c(0, 1000))
lines(run2[[1]], col="red")
lines(run3[[1]], col="blue")
lines(run4[[1]], col="green")

#Theta
plot(run1[[2]], type="l", col="black", xlab="Iterations", ylab="Theta", ylim=c(0, 0.00001), xlim=c(0, 1000))
lines(run2[[2]], col="red")
lines(run3[[2]], col="blue")
lines(run4[[2]], col="green")

#Impact of the priors on the posterior_tau
sample_tau = run1[[1]][50:length(run4[[1]])];
summary(sample_tau)  
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run2[[1]][50:length(run4[[1]])];
summary(sample_tau)  
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run3[[1]][50:length(run4[[1]])];
summary(sample_tau)  
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

sample_tau = run4[[1]][50:length(run4[[1]])];
summary(sample_tau)  
hist(sample_tau, main="posterior", xlab="Tau", freq=FALSE, n=20)
mean_tau <- mean(sample_tau)
mean_tau

#Impact of the priors on the posterior_theta
sample_theta = run1[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run2[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run3[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

sample_theta = run4[[2]][800:length(run4[[2]])];
summary(sample_theta)  #posterior mean for 1st MCMC, dropping 500 values as burn-in
hist(sample_theta, main="posterior", xlab="Theta", freq=FALSE, n=20)
mean_theta <- mean(sample_theta)
mean_theta

#L

L <- c(1,5,10,50,100,150,175,200,250,300,400,500,600,700,800,900,1000,5000,10000)
theta <- c(0.001089818, 0.001120797, 0.001061077,0.001049477,0.001074818,0.0008421703,0.0004156964,0.0002774716, 0.0001255517,4.730733e-05, 6.380365e-06, 1.846478e-05, 4.454329e-07, 3.982985e-06, 6.787579e-06, 2.300075e-06, 5.512387e-07, 2.49454e-06, 8.890647e-07)

plot(L, theta, main = "Posterior",
     xlab = "L", ylab = "Theta", ylim=c(0, 0.0015), xlim=c(0, 1000),
     pch = 19, frame = FALSE)



