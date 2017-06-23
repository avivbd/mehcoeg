
graphics.off()  #close plots
rm(list = ls()) #clear global env
cat("\014") #clear console



#define our "unkown" function
unkownfun1 = function(x){
  y = dgamma(x,2,1)
   if(x > 2){
      y = 0
   }
  return(y)
  #dnorm(x,mean = 3, sd = 2)
}


# the MCMC function. The arguments being passed to it are:
# n the number of iterations, defualt = 1000.
# the min/max of the uniform distribution from which proposed steps are drawn, default = 0.5
# the initial guess, defualt = 0
# a call to the unkown function that we're trying to recreate
# a call to the proposal distribution (the distribution from which proposals are drawn)

propdist1 <- function(x, mu = 0 , eps = 1){
  #CDF = pnorm(x) #draw from some dist
  CDF = pt(x, 5) #draw from some dist
  #   #dense = dexp(x)  #draw from some dist
  #   dense = runif(1, -eps, eps) #draw from some dist
  return(CDF)
}

randraw1 <- function(x){
  #innov = rnorm(1)#
  innov = rt(1,5)#
  can = x + innov # candidate value to evaluate function at
  return(can)
}

metrop1 <- function(n = 1000, eps = 1, x0 = 0, unknfun, propdist) 
{
  
  x <- x0 # initial guess
  oldlik = unknfun(x) # first function evaluation
  
  vec = vector("numeric", n) # preallocate memory for unkown parameter values
  vec[1] = x # put in the initial guess
  
  liklihood <- vector("numeric",n)
  liklihood[1] = oldlik
  
  for (i in 2:n){
    #repeat {
      #xs=x+rnorm(1)
      can = randraw1(x)
     # if (can>0)
      #  break
  #  }
    
    #calculate proposal correction factor
    c = propdist(x)/propdist(can)
    # c = propdist(x, can, eps)/propdist(can, x, eps)
    #print(can)
    lik = unknfun(can) # value of unkown function at candidate location
    a = lik/oldlik#*c # ratio of the function values
    if(is.na(a)){a = 0}
    #print(a)
    u = runif(1) # random number on [0 1]
    if (u < a){ # test if liklihood ratio is larger than random number. 
      # If candidate location is better than old location
      # the proposal will always be accepted. If worse 
      # will be accepted sometimes.
      x = can
      oldlik = lik
      #points(x,lik) #draw points on likelihood function. Will slow down code. comment out if you don't want it. 
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
  #since we actually know it we can plot it in advance and see how well MCMC recreates it
  
  plot(ts(mcmc.out),col=2)
  hist(mcmc.out,30,col=3,freq = F)
  x = seq(-6,6,length=12*4)
  f = unkownfun1(x)
  lines(x,f,type = "l")
  
  # qqnorm(mcmc.out,col=4)
  # abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}

mcmc_out <- metrop1(n = 10000, x0 = -1, eps= 1, unknfun = unkownfun1, propdist = propdist1)
plot_mcmc(mcmc_out)
