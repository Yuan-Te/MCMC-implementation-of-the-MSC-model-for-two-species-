### MCMC program in R for estimating JC69 distance
### R program implementing an MCMC algorithm for distance estimation under JC69 model ### with exponential prior. This is for the example on page 163 and figure 5.7 of ### Yang (2006 Computational Molecular Evolution, OUP).
### Code written by Konstantinos Angelis <aggelisk@gmail.com>, edited by Ziheng Yang ### 3 February 2013.
### Code is free. Please copy and edit as you like.
logpriorlikelihood = function(t, x, n, mu){
  # logarithm of likelihood and prior for JC69 distance.  We calculate the log
  # to avoid underflows and overflows when t, x, or n are large.
  # x = number of observed nucleotide differences between two sequences
  # n = sequence length
  # t = distance
  # mu = mean of exponential prior for t
  p = 3/4 - (3/4)*exp(-(4*t)/3) # p也等於 x/n
  lnP = x*log(p) + (n-x)*log(1-p) -log(mu) - t/mu
  return(lnP) #lnP is log(posterior)
}


# The term "-log(mu)" is not necessary
mcmc = function(t, N, x, n, window, mu){
  # t(Initial) = initial value of distance.  t & tnew are the distance parameter
  # N = number of samples including initial state
  # x = number of observed differences between two sequences
  # n = sequence length
  # window = window size
  # mu = prior mean of distance
  sample = numeric(N+1)
  sample[1] = t
  accept = 0   #counts the number of acceptances
  lnp = logpriorlikelihood(t, x, n, mu);
  
  for (i in 1:N) {
    tnew = t - window/2 + window*runif(1);
    if (tnew < 0) tnew = -tnew;  # reflection at boundary t=0
    lnpnew = logpriorlikelihood(tnew, x, n, mu);
    logratio = lnpnew - lnp;
    
    if(logratio>=0 || runif(1)<exp(logratio)) {  # accept tnew
      t = tnew;
      lnp = lnpnew;
      accept = accept+1
    }
    sample[i+1] = t
  }
  accept = accept/(N - 1)
  return(list(sample, accept))
}  #end of function


#-------------------------------------------------------------------------------------------

jc_run1 = mcmc(t=0.000001, N=20000, x=90, n=948, window=0.1, mu=0.2)  #initial value = 0.01
jc_run2 = mcmc(t=0.05,  N=20000, x=90, n=948, window=0.1, mu=0.2)  #initial value = 0.5
jc_run3 = mcmc(t=1,    N=20000, x=90, n=948, window=0.1, mu=0.2)  #initial value = 1

#accept rate for the 3 MCMC runs
c(run1[[2]], run2[[2]], run3[[2]])

plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="Distance", ylim=c(0, 1), xlim=c(0, 200))
lines(run2[[1]], col="red")
lines(run3[[1]], col="blue")

sample1 = jc_run1[[1]][101:20000];
summary(sample1)  #posterior mean for 1st MCMC, dropping 100 values as burn-in
hist(sample1, main="posterior", xlab="distance", freq=FALSE, n=40)

#exercise:
#run a chain with tInitial=1, N=1e5, window=0.01
#run a chain with tInitial=1, N=1e5, window=0.1
#run a chain with tInitial=1, N=1e5, window=10
#plot the three runs in one graph
#which window sizes are good and which are bad?


#auto-correlation function & efficiency
#
sample1.acf = acf(sample1)   # The correlation is lost after a lag of k=10
sample1.acf$acf              # the actual acf values are here.


# If the autocorrelation is high, the chain will be inefficient, i.e. we will have to run a very 
# long chain to obtain a good approximation to the posterior.
# The efficiency of a chain is defined as:
# eff = 1 / (1 + 2(r1 + r2 + r3 + ...))
# where r_k is the correlation for lag k.

eff = function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))
eff(sample1.acf)  # ~ 25% What does that mean?

x
n
