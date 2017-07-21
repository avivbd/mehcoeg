library(mvtnorm)

#The R functions for the particle/bootstrap filter (PF) and the
#backward-simulation particle smoother can be downloaded from: 
# source("http://www.datall-analyse.nl/R/PF.R")
source("/Users/avivb/Google Drive/Research/StanfordPostdoc/KalmanFiltering/code/R/PF.R")
#Take a look at the function PF (=the particle/bootstrap filter) and
#PFsmooth (=smoother), and notice that at the beginning of the scripts
#you will find a description of the functions' arguments.
# PF
# PFsmooth



# EXAMPLE 1 ----

#The following example can be found in Crassidis and Junkins's
#2012 book "Optimal estimation of dynamic systems".

#In this example we will estimate the state of the following
#nonlinear discrete-time system:
#x(k+1) = x(k)/2 + 25*x(k)/(1+x(k)^2) + 8*cos(1.2*k) + W
#where x(k+1) is the state at time step k+1, x(k) the state at time
#step k, and W the process noise (variance). We will set W to 10 in
#this example.

#The measurement equation is:
#y(k) = x(k)^2/20 + V
#where V is the measurement noise (variance). We will set V to 1.



##TRUE STATES AND MEASUREMENTS

#sample time
dT <- 1

#observation times
t <- seq(1, 75, dT)

#process noise (standard deviation)
Wk <- sqrt(10)

#measurement noise (standard deviation)
Vk <- sqrt(1)

#state at t=0
X0 <- 0 + rnorm(1)*sqrt(5)

#simulate true state for x at the observation times
Xk <- matrix(NA, nrow=length(t)+1, ncol=1)
Yk <- matrix(NA, nrow=length(t), ncol=1)
Xk[1] <- X0

for (i in 2:(length(t)+1)) {
  Xk[i] <- Xk[i-1]/2 + (25*Xk[i-1])/(1+Xk[i-1]^2) + 8*cos(1.2*(i-1)) + rnorm(1)*Wk
  Yk[i-1] <- (Xk[i]^2)/20 + rnorm(1)*Vk}

#plot simulated states
plot(c(0,t), Xk, type="o", xlab="Time, k", ylab="x(k)")

#plot measurements
plot(t, Yk, type="o", xlab="Time, k", ylab="y(k)")

#store measurement data
dataEx1 <- Yk



##PARTICLE/BOOTSTRAP FILTER (PF)

#Dynamic model:
#specifying 1 state, namely [x].
#Note that you may change the initial state estimate (at t=0) below for x
#and see how it influences the behavior of the particle filter.
ex1 <- list(m0=0, #initial state estimate
            #error covariance of the initial state estimate:
            #this vector reflects the uncertainty in our initial state estimate
            #you may change the value in this vector and see how it influences
            #the behavior of the particle filter
            C0=5,
            #measurement noise
            V=Vk^2,
            #process noise
            W=Wk^2)

#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  x/2 + (25*x)/(1+x^2) + 8*cos(1.2*(k))}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  (x^2)/20}


##Compute the filtered (a posteriori) state estimates with the PF
pf1 <- PF(y=dataEx1, mod=ex1, N=500, resampling="strat",
          GGfunction=GGfunction, FFfunction=FFfunction)



##Perform backward-simulation particle smoothing.
#Simulate 20 trajectories:
pfs <- replicate(20, PFsmooth(pf1))

#Plot the smoother results for the 20 simulated trajectories
pfsm <- do.call(cbind, pfs)

plot(c(0, t), Xk, type="o", col="black", lwd=2, ylab="x(k)", xlab="Time, k")
matlines(c(0, t), pfsm, lty=2, col=gray(level=.5, alpha=.8), lwd=1)
lines(c(0, t), pf1$m[,1], lty=2, col="blue", lwd=2)
legend("topleft", col=c("black", "blue", gray(level=.5)),
       lty=c(1, 2, 2), pch=c(1, NA, NA), lwd=c(2, 2, 1),
       legend=c("true state", "filtered PF", "smoothed PF"),
       bty="n", y.intersp=1.2, cex=.7)


#Average of 20 simulated trajectories
plot(c(0, t), Xk, type="o", col="black", lwd=2, ylab="x(k)", xlab="Time, k")
points(c(0, t), apply(pfsm, 1, mean), lty=2, col="red", lwd=1)
lines(c(0, t), pf1$m[,1], lty=2, col="blue", lwd=1)
legend("topleft", col=c("black", "blue", "red"),
       lty=c(1, 2, NA), pch=c(1, NA, 1), lwd=c(2, 1, 1),
       legend=c("true state", "filtered PF", "smoothed PF"),
       bty="n", y.intersp=1.2, cex=.7)