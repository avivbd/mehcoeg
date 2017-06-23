##  #parallel tempering and Metropolis MCMC


#demo 1
U=function(gam,x)
{
  gam*(x*x-1)*(x*x-1)
}

curried=function(gam)
{
  message(paste("Returning a function for gamma =",gam))
  function(x) U(gam,x)
}



op=par(mfrow=c(2,1))
curve(curried(20)(x),-2,2,main="Potential function, U(x)") #cost function
curve(exp(-curried(20)(x)),-2,2,main="Unnormalised density function, exp(-U(x))") #liklihood function
par(op)

#demo 2 unlinked chains running sequentially
#Metropolis Algorithm
chain=function(target,tune=0.1,init=1)
{
  x=init
  xvec=numeric(iters)
  for (i in 1:iters) {
    can=x+rnorm(1,0,tune)
    logA=target(x)-target(can)
    if (log(runif(1))<logA)
      x=can
    xvec[i]=x
  }
  xvec
}

temps=2^(0:3)
iters=1e5

mat=sapply(lapply(temps,curried),chain)
colnames(mat)=paste("gamma=",temps,sep="")

require(smfsb)
mcmcSummary(mat,rows=length(temps))


#demo 3 parallel chains

chains=function(pot=U, tune=0.1, init=1)
{
  x=rep(init,length(temps))
  xmat=matrix(0,iters,length(temps))
  for (i in 1:iters) {
    can=x+rnorm(length(temps),0,tune)
    logA=unlist(Map(pot,temps,x))-unlist(Map(pot,temps,can))
    accept=(log(runif(length(temps)))<logA)
    x[accept]=can[accept]
    xmat[i,]=x
  }
  colnames(xmat)=paste("gamma=",temps,sep="")
  xmat
}

mcmcSummary(chains(),rows=length(temps))

#demo 4 coupled chains with swapping

chains=function(pot=U, tune=0.1, init=1)
{
  x=rep(init,length(temps))
  xmat=matrix(0,iters,length(temps))
  for (i in 1:iters) {
    can=x+rnorm(length(temps),0,tune)
    logA=unlist(Map(pot,temps,x))-unlist(Map(pot,temps,can))
    accept=(log(runif(length(temps)))<logA)
    x[accept]=can[accept]
    # now the coupling update
    swap=sample(1:length(temps),2)
    
    logA=pot(temps[swap[1]],x[swap[1]])+pot(temps[swap[2]],x[swap[2]])-
      pot(temps[swap[1]],x[swap[2]])-pot(temps[swap[2]],x[swap[1]])
    #careful this is confusing: uses -U(x) so the swapped is on the bottom
    
    if (log(runif(1))<logA)
      x[swap]=rev(x[swap])
    # end of the coupling update
    xmat[i,]=x
  }
  colnames(xmat)=paste("gamma=",temps,sep="")
  xmat
}


mcmcSummary(chains(),rows=length(temps))









# 
# chain=function(target,tune=0.1,init=1)
# {
#   x=init
#   xvec=numeric(iters)
#   for (i in 1:iters) {
#     can=x+rnorm(1,0,tune)
#     logA=target(x)-target(can)
#     if (log(runif(1))<logA)
#       x=can
#     xvec[i]=x
#   }
#   xvec
# }
# 
# 
# temps=2^(0:3)
# iters=1e5
# 
# mat=sapply(lapply(temps,curried),chain)
# colnames(mat)=paste("gamma=",temps,sep="")
# 
# require(smfsb)
# mcmcSummary(mat,rows=length(temps))
# 
# 
# 
# chains2=function(pot=U, tune=0.1, init=1)
# {
#   x=rep(init,length(temps))
#   xmat=matrix(0,iters,length(temps))
#   for (i in 1:iters) {
#     can=x+rnorm(length(temps),0,tune)
#     logA=unlist(Map(pot,temps,x))-unlist(Map(pot,temps,can))
#     accept=(log(runif(length(temps)))<logA)
#     x[accept]=can[accept]
#     xmat[i,]=x
#   }
#   colnames(xmat)=paste("gamma=",temps,sep="")
#   xmat
# }
# 
# mcmcSummary(chains2(),rows=length(temps))
