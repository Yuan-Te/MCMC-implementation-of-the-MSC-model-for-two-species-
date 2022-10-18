lnPf <- function(p,q,r){
  #Posterior = Prior x likelihood
  #ln(posterior) = ln(prior)+ln(likelihood)
  ln_posterior = log(2) + 44*log(p^2+2*p*r) + 27*log(q^2+2*p*r) + 4*log(2*p*r) + 88*log(r^2)
  return(ln_posterior)
}

ln_posterior <- lnPf(p=0.3, q=0.3, r=0.4) #Test
posterior <- exp(ln_posterior)

#Reflection
reflect <- function(x,a,b){
  #This returns the value when x is reflected into the range(a,b)
  
  side = 0; #side = 0(left) or 1(right).
  e = 0; #excess
  
  if(x<a){
    e = a-x; side = 0;
  }else if(x>b){
    e = x-b; side = 1;
  }
  
  if(e != 0){
    n = trunc(e/(b-a));
    
    if(n%%2 != 0){
      side = 1 - side;
    }
    
    e = e - n*(b-a);
    if(side == 1){
      x = b-e
    }else {
      x = a+e
    } 
  }
  
  return(x)
} 

test <- reflect(-1.8,0,1) #Test

abo = function(N, p, q, r, w){
  # tInitial = initial value of distance. t & tnew are the distance parameter
  # N = number of samples
  # x = number of observed differences between two sequences
  # n = sequence length
  # window = window size
  # mu = prior mean of distance
  
  if(p<0 || q<0 || r<0 || abs(1-p-q-r)>1e-6) {
    print("bad initial values");
    return (0);
  }
  
  lnP = lnPf(p, q, r)
  sample_p = numeric(N+1)
  sample_q = numeric(N+1)
  sample_r = numeric(N+1)
  sample_p[1] = p
  sample_q[1] = q
  sample_r[1] = r
  
  accept = 0
  
  for(i in 1:N){
    u = runif(1) #创建均匀分布的随机偏差, 大U的意思
    if(u<1/3){
      s = p+q   # change p and q with s = p + q fixed
      pnew = p + w*(runif(1)-0.5)
      pnew = reflect(pnew, 0, s)
      qnew = s - pnew
      rnew = r
      
    }else if(u<2/3){
      s = q+r
      qnew = q + w*(runif(1)-0.5)
      qnew = reflect(qnew, 0, s)
      rnew = s - qnew
      pnew = p
      
    }else{
      s = r+p
      rnew = r + w*(runif(1)-0.5)
      rnew = reflect(rnew, 0, s)
      pnew = s - rnew
      qnew = q
    }
    
    lnPnew = lnPf(pnew, qnew, rnew)
    lnalpha = lnPnew - lnP
    if(lnalpha>0 || exp(lnalpha) > runif(1)){  
      p = pnew; q = qnew; r = rnew; lnP = lnPnew
      accept = accept+1
    }
    #If π(θ*) > π(θ), accept the proposal; otherwise accept it with probability π(θ*) / π(θ).
    #Generating a random number r from U(0, 1), if r < α, accept θ*, or otherwise reject θ*.
    #This way the proposal is accepted with probability α.
    
    sample_p[i+1] = p
    sample_q[i+1] = q
    sample_r[i+1] = r
  }
  
  accept = accept/N
  return(list(sample_p, sample_q, sample_r, accept))
}

test <- abo(N=2, p=0.3, q=0.3, r=0.4, w=0.005) #Test


run1 = abo(N=20, p=0.3, q=0.3, r=0.4, w=0.005)
run1[[4]] #accept rate
plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="p",ylim=c(0,1)) #sample_p
lines(run1[[2]], col="red") #sample_q
lines(run1[[3]], col="blue") #sample_r

run1 = abo(N=200, p=0.3, q=0.3, r=0.4, w=0.1)
run1[[4]] #accept rate
plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="p",ylim=c(0,1)) #sample_p
lines(run1[[2]], col="red") #sample_q
lines(run1[[3]], col="blue") #sample_r

run1 = abo(N=10000, p=0.3, q=0.3, r=0.4, w=0.08)
run1[[4]] #accept rate
plot(run1[[1]], type="l", col="black", xlab="Iterations", ylab="p",ylim=c(0,1)) #sample_p
lines(run1[[2]], col="red") #sample_q
lines(run1[[3]], col="blue") #sample_r
    
summary(run1[[1]])
summary(run1[[2]])
summary(run1[[3]])

plot(run1[[1]], run1[[2]]) #不懂

hist(run1[[1]], breaks=50)
hist(run1[[2]], breaks=50)
hist(run1[[3]], breaks=50)

#https://easystats.github.io/bayestestR/articles/credible_interval.html
library(bayestestR)
library(dplyr)
library(ggplot2)
ci_hdi <- ci(run1[[1]], method = "HDI")
ci_eti <- ci(run1[[1]], method = "ETI")
run1[[1]] %>% 
  estimate_density(extend=TRUE) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_area(fill = "orange") +
  theme_classic() +
  # HDI in blue
  geom_vline(xintercept = ci_hdi$CI_low, color = "royalblue", size = 3) +
  geom_vline(xintercept = ci_hdi$CI_high, color = "royalblue", size = 3) +
  # Quantile in red
  geom_vline(xintercept = ci_eti$CI_low, color = "red", size = 1) +
  geom_vline(xintercept = ci_eti$CI_high, color = "red", size = 1)

ci_hdi <- ci(run1[[2]], method = "HDI")
ci_eti <- ci(run1[[2]], method = "ETI")
run1[[2]] %>% 
  estimate_density(extend=TRUE) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_area(fill = "orange") +
  theme_classic() +
  # HDI in blue
  geom_vline(xintercept = ci_hdi$CI_low, color = "royalblue", size = 3) +
  geom_vline(xintercept = ci_hdi$CI_high, color = "royalblue", size = 3) +
  # Quantile in red
  geom_vline(xintercept = ci_eti$CI_low, color = "red", size = 1) +
  geom_vline(xintercept = ci_eti$CI_high, color = "red", size = 1)



