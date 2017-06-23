
graphics.off()  #close plots
rm(list = ls()) #clear global env
cat("\014") #clear console


#make a time-series 
tt = seq(from = 0, to = 10, length.out = 1000)
slope = 2
intercept = -1
#model function whose parameters we are interested in 
modelfun = function(slope, intercept){
  Y = slope%*%tt + intercept
  return(Y)
}  
#noisy data
Data = modelfun(slope, intercept)  + rnorm(length(tt)) 

#cost function evaluates the root mean square error normalized by the number of datapoints
costfun = function(unkownparams){
  #model output should be horizontal in time
  #and different rows should be different iterations
  sizeY = dim(modelfun(unkownparams(1),unkownparams(2)))
  NRMSE = vector("numeric",length = sizeY[1])
   for (i in 1:sizeY[1]){
    NRMSE[i] = sqrt(  sum((modelfun(unkownparams[i]) - Data)^2)/sizeY[2] )
     #/( max(TS - x) - min(TS - x) )
   }
   return(NRMSE)
}

#liklihood function returns the negative exponent of the cost function
likefun = function(x){
  likefun = exp(-costfun(x))
  return(likefun)
}


#for testing
# par(mfrow = c(3,1))
# plot(tt,Data,ylab = "data")
# plot(tt, modelfun(slope, intercept),ylab = "model")
# unkownmu = -10:15
# plot(unkownmu, likefun(unkownmu),xlab = "unkown param",ylab = "neg log cost fun")


#metropolis algorithm
#metrop3 <- function(n = 10000,eps = .1, x0 = rnorm(2) ){
  eps = 0.1
  n = 10000
  x = c(1,1)
  vec <- vector("numeric",n)
  lik <- vector("numeric",n)
  x = x0
  OldLoglik = log(likefun(x))
  for (i in 2:n){
    #innov = runif(1,-eps,eps)
    innov = rnorm(1,mean = 0,sd = eps)
    can = x + innov
    CanLoglik = log(likefun(can))
    loga = CanLoglik - OldLoglik
    logu = log(runif(1))
    if (logu < loga){
      x <- can
      OldLoglik <- CanLoglik
    }
    vec[i] = x
    lik[i] = OldLoglik
  }
  #return(rbind(vec,lik))
#}




N_MCMC_Iterations = 20000
mcmc_out <- metrop3(n = N_MCMC_Iterations)
theta = mcmc_out[1,]
liklihood = mcmc_out[2,] 
#plotting and statistics
op = par(mfrow= c(2,3))
plot(1:N_MCMC_Iterations,theta,type = "l", col = 1,main = "Markov Chain")
burntime = ceiling(0.1*N_MCMC_Iterations)
theta_no_burn <- theta[-1:-burntime]
hist(theta_no_burn,30,col=3,freq = F,xlab = "value of unkown param")
#acf(mcmc_out,col=2,lag.max=100)
qqnorm(theta_no_burn,col=4)
qqline(theta_no_burn,col=2)
#plot the data and the model
plot(tt,Data,xlab = "time", ylab = "data", main = "Plot of model and data")
unkown_density = density(theta_no_burn)
MostLikelyParamValue = unkown_density$x[which.max(unkown_density$y)]
lines(tt, modelfun(MostLikelyParamValue),col = 2)

  
  plot(theta, exp(liklihood),main = "likelihood function" )

