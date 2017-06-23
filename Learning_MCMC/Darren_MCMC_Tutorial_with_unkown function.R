
graphics.off()  #close plots
rm(list = ls()) #clear global env
cat("\014") #clear console



#define our "unkown" function
unkownfun1 = function(x){
  dnorm(x)
}

#since we actually know it we can plot it in advance and see how well MCMC recreates it
x11()
x = seq(-6,6,length=12*4)
f = unkownfun1(x)
plot(x,f,type = "l")

# the MCMC function. The arguments being passed to it are:
# n the number of iterations, defualt = 1000.
# the min/max of the uniform distribution from which proposed steps are drawn, default = 0.5
# the initial guess, defualt = 0
# a call to the unkown function that we're trying to recreate
metrop1 <- function(n = 1000, eps = 0.5, x0 = 0, myfun) 
{
  
  x <- x0 # initial guess
  oldlik = myfun(x) # first function evaluation
  
  vec = vector("numeric", n) # preallocate memory for unkown parameter values
  vec[1] = x # put in the initial guess
  
  liklihood <- vector("numeric",n)
  liklihood[1] = oldlik
  
  for (i in 2:n){
    innov = runif(1,-eps,eps) #draw from random uniform dist
    can = x + innov # candidate value to evaluate function at
    lik = myfun(can) # value of unkown function at candidate location
    a = lik/oldlik # ratio of the function values
    u = runif(1) # random number on [0 1]
    if (u < a){ # test if liklihood ratio is larger than random number. 
      # If candidate location is better than old locattion
      # the proposal will always be accepted. If worse 
      # will be accepted sometimes.
      x = can
      oldlik = lik
      points(x,lik) #draw points on likelihood function. Will slow down code. comment out if you don't want it. 
    }
    vec[i] = x #save values
    liklihood[i] = lik #save liklihood
  }
  return(vec)
  
}


#plotting function
plot_mcmc<-function(mcmc.out)
{
  op=par(mfrow=c(2,2))
  plot(ts(mcmc.out),col=2)
  hist(mcmc.out,30,col=3,freq = F)
  qqnorm(mcmc.out,col=4)
  abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}

mcmc_out <- metrop1(n = 5000, eps = 1, x0 = 5, myfun = unkownfun1)
plot_mcmc(mcmc_out)
