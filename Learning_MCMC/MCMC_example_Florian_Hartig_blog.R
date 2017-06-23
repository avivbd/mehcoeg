# MCMC from Theoretical ecology blog Florian  Hartig

graphics.off()  #close plots
rm(list = ls()) #clear global env
cat("\014") #clear console

#define mystery function which underlies the noisy data
trueA = 5
trueB = 0
trueSD = 10
sampleSize = 31

x = (-(sampleSize-1)/2):((sampleSize-1)/2)
y = trueA*x + trueB + rnorm(n = sampleSize, mean = 0, sd = trueSD)
plot(x,y, main = 'test data')

#Log likelihood function. Gives -125.1104 for True values
likelihood <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  pred = a*x + b
  
  singlelikelihoods = dnorm(y,mean = pred, sd = sd, log = T)
  sumll = sum(singlelikelihoods)
  
  #alternative log likelihood function. Gives -162 for real values
  #log(exp(-sum((y-pred)^2)/(length(y)))) = -sum((y-pred)^2)/(length(y))
  #sumll = -sum((y-pred)^2)/(length(y))
  
  return(sumll)
}


slopevalues <- function(x){
  return(likelihood(c(x,trueB,trueSD)))
}
xx = seq(3,7,by = 0.05)
slopelikelihoods = lapply(xx,slopevalues)
plot(xx, slopelikelihoods)

#defining priors
prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  aprior = dunif(a, min = 0, max = 10, log=T)
  bprior = dnorm(b, sd = 5, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+sdprior)
}

#posterior is the product of the prior and likelihood or their sums in log space
posterior <- function(param){
  return( likelihood(param) + prior(param) )
}
  

proposalfunction <- function(param){
  y = rnorm(3, mean = param, sd = c(0.1, 0.5, 0.3))
  return(y)
}
  
run_metropolis_mcmc <- function(startvalue,iterations){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}
  
startvalue = c(4,0,10)
chain = run_metropolis_mcmc(startvalue, 100000)
  
burnIn = 5000
rejection = 1 - mean(duplicated(chain[-(1:burnIn),]))

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = trueA, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = trueB, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSD, col="red" )
plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = trueA, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = trueB, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSD, col="red" )

# for comparison:
summary(lm(y~x))

#look at the joints
pairs(data.frame(chain))

