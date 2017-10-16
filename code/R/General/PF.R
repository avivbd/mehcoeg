PF <- function (y, mod, GGfunction, FFfunction, N,
                Nthreshold=N+1, resampling="strat",
                roughening=FALSE, Grough=.2,
                MCparticles=TRUE, logLik=FALSE, simplify=FALSE) {
  
  #Originally written by Stephan Gelissen. 
  #http://blogs2.datall-analyse.nl/2016/04/11/rcode_particle_bootstrap_filter/    
  #Modified by Aviv Bachan.     
    
  #This function implements the particle filter (PF) as described
  #by Crassidis and Junkins in their 2012 book "Optimal estimation of
  #dynamic systems". This PF is also referred to as the bootstrap filter.
  
  
  #Note that this PF assumes (time invariant) additive Gaussian noise.
  
  #Last update PF function: 2016/06/16
  #Modified by AB 07/17
  
  #Arguments:
  #-y is the data. The data can be a vector (one-dimensional measurements /
  # univariate time series) or a matrix (multidimensional measurements / 
  # multivariate time series). In the matrix case, each of the columns contains
  # the measurements on a single dimension (or, similarly, a single [univariate]
  # time series).
  # Missing values (NA's) are allowed.
  #-mod is a list with the components m0 (initial state estimates),
  # C0 (error covariances of the initial state estimates), V (measurement noise),
  # and W (process noise).
  #-GGfunction is a function with arguments x and k. The GGfunction specifies the
  # state transition function. The transition function is usually denoted as
  # f(x[k], u[k]), where x are the states estimates and u the control input at
  # time step k.
  # See details below for how to handle control input in the GGfunction.
  #-FFfunction is a function with arguments x and k. The FFfunction specifies the
  # observation/measurement function. The measurement function is usually denoted
  # as h(x[k]).
  #-N is the number of particles.
  #-Nthreshold is a number specifiying the resampling threshold. That is, resampling
  # takes place when the effective sample size falls below the threshold.
  #-resampling is the resampling method to be applied: either "mult"
  # (=multinomial resampling), "strat" (=stratified resampling), or
  # "sys" (=systematic resampling).
  #-roughening: if TRUE, then perform roughening.
  #-Grough is a number specifying the roughening tuning parameter.
  #-MCparticles: if TRUE, then include the state estimates (=xpt), and
  # the importance weights (=wt) of each particle at each time step in the output.
  #-logLik: if TRUE, then return the negative loglikelihood of the specified model.
  #-simplify: if TRUE, then do not include the data (=y) in the output.
  
  #Details:
  #It is possible to include control input in the GGfunction.
  #As an example, let us assume that there are two states (x1 and x2) and that
  #we have 70 time steps at which we want to obtain state estimates. At each of
  #these 70 time steps there will be (external) control input.
  #
  #We may specify an arbitrary GGfunction with control input as follows:
  #GGfunction <- function (x, k){
  #x1 <- x[1]; x2 <- x[2]
  #c(2*x1 + 1*x2 - 4,
  #  3*x2 - 1*x2 - 2)}
  #where, at each time step, the control input for the first state equation
  #is -4, and the control input for the second state equation -2.
  #
  #Alternatively, we may store the control input for the two state equations
  #at each of the 70 time steps in a 70x2 matrix.
  #Subsequently, we specify the GGfunction as:
  #GGfunction <- function (x, k){
  #x1 <- x[1]; x2 <- x[2]
  #c(2*x1 + 1*x2 - u[k,1],
  #  3*x2 - 1*x2 - u[k,2])}
  #where we have stored the control input at each of the 70 time steps
  #in a matrix called u.
  #
  #Note that the latter specification of the GGfunction is useful when the control
  #input is known to vary over time.
  
  
    
  mod1 <- mod
  y <- as.matrix(y)
  ym <- ncol(y)
  yAttr <- attributes(y)
  p <- length(mod$m0)
  
  if (!is.null(mod$FF) | !is.null(mod$GG))
    warning ("FF or GG matrix will not be used in the PF")
  
  if (!is.null(mod$JFF) | !is.null(mod$JGG))
    warning ("Time varying FF or GG matrix will not be used in the PF")
  
  if (!is.null(mod$JW) | !is.null(mod$JV))
    warning ("Time varying V or W matrix will not be used in the PF")
  
  m <- rbind(mod$m0, matrix(0, nrow = nrow(y), ncol = length(mod$m0)))
  a <- matrix(0, nrow = nrow(y), ncol = length(mod$m0))
  f <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  C <- vector(1 + nrow(y), mode = "list")
  R <- vector(nrow(y), mode = "list")
  if (MCparticles) {
    #importance weights before resampling at each time step
    wt <- matrix(0, nrow = nrow(y) + 1, ncol = N)
    #particles at each time step
    xpt <- vector(1 + nrow(y), mode = "list")}
  ll <- 0
  invV <- solve(mod$V)
  
  
  #generate particles at t=0
  xp <- rmvnorm(n=N, mean=m[1,], sigma=as.matrix(mod$C0))
  
  # print(xp)
  # browser()
  
  #importance weights at t=0
  w <- rep(1/N, N)
  
  if (MCparticles) {
    wt[1, ] <- w
    xpt[[1]] <- xp}
  
  C[[1]] <- mod$C0
  
  for (i in seq(length = nrow(y))) {
    
    if (!any(whereNA <- is.na(y[i, ]))) {
      
      ##time update
      wks <- rmvnorm(n=N, sigma=as.matrix(mod$W))
      #a priori state estimate
      model_predict = apply(xp, 1, GGfunction, k=i)
      xp <- t(matrix(model_predict, nrow=p)) + wks
      
      
      a[i, ] <- crossprod(w, xp)
      #covariance of a priori state estimate
      xc <- t(t(xp)-a[i, ])
      R[[i]] <- crossprod(sqrt(w)*xc)
      
      ##measurement update
      vks <- rmvnorm(n=N, sigma=as.matrix(mod$V))
      #predicted measurement
      obs_predict=apply(xp, 1, FFfunction, k=i)
      yp <- t(matrix(obs_predict, nrow=ym)) + vks
      
      # print(c("mod_pred: ",colMeans((t(model_predict)))))
      # print(c("mod_pred+W: ",colMeans(((xp)))))
      # print(c("mean W %: " , colMeans(abs(wks))/mod$m0*100 ))
      # print(c("obs_pred: ",colMeans((t(obs_predict)))))
      # print(c("obs_pred+V: ",colMeans((yp))))
      # print(c("mean V %: ", colMeans(abs(vks))/y[i,]*100))
      # print(c("obs: ", y[i,]))
      
      
      f[i, ] <- crossprod(w, yp)
      #covariance of predicted measurement
      yc <- t(t(yp)-f[i, ])
      Qy <- crossprod(sqrt(w)*yc)
      
      ##a posteriori estimates
      
      #importance weights update
      xs <- matrix(y[i, ] - t(yp), nrow=N, byrow=TRUE)
      w <- w*apply(xs, 1, function (x) exp(-.5*(tcrossprod(crossprod((x), invV), t(x)))))
      if (all(w == 0)) next() #Aviv's hack. Won't die if all the weights are zero. 
      w <- w/sum(w)
      
      
      # summary(xp) %>% print
      # summary(yp) %>% print
      # browser()
      
      #a posteriori state estimate
      m[i + 1, ] <- crossprod(w, xp)
      #a posteriori error covariance
      C[[i + 1]] <- crossprod(sqrt(w)*xc)
      
      if (MCparticles) {
        xpt[[i + 1]] <- xp
        wt[i + 1, ] <- w}
      
      ##resampling
      Neff <- 1 / crossprod(w)
      
      if (Neff < Nthreshold) {
        
        if (resampling=="mult") {
          xp <- as.matrix(resampleM(n=N, wk=w, xk=xp))
          w <- rep(1/N, N)}
        else if (resampling=="sys") {
          xp <- as.matrix(resampleSy(n=N, wk=w, xk=xp))
          w <- rep(1/N, N)}
        else if (resampling=="strat") {
          xp <- as.matrix(resampleSr(n=N, wk=w, xk=xp))
          w <- rep(1/N, N)}
      }
      
      ##roughening
      if (roughening) {
        El <- apply(xp, 2, function (x) diff(range(x)))
        sigma.l <- Grough*El*N^(-1/length(mod$m0))
        if (length(El)==1) Jk <- sigma.l^2 else
          Jk <- diag(sigma.l^2)
        ck <- rmvnorm(n=N, sigma=as.matrix(Jk))
        xp <- xp + ck}
      
      
      # summary(xp) %>% print
      # summary(yp) %>% print
      # browser()
      
      ##compute log-likelihood
      if (logLik) {
        e <- as.matrix(y[i, ] - f[i,])
        ll <- ll + ym*log(2*pi) + sum(log(eigen(Qy)$values)) + 
          crossprod(e, tcrossprod(solve(Qy, tol=1e-30), t(e)))}
    }
    
    else {
      if (all(whereNA)) {
        
        ##time update
        wks <- rmvnorm(n=N, sigma=as.matrix(mod$W))
        #a priori state estimate
        xp <- t(matrix(apply(xp, 1, GGfunction, k=i), nrow=p)) + wks
        a[i, ] <- crossprod(w, xp)
        #covariance of a priori state estimate
        R[[i]] <- crossprod(sqrt(w)*t(t(xp)-a[i, ]))
        
        ##measurement update
        vks <- rmvnorm(n=N, sigma=as.matrix(mod$V))
        #predicted measurement
        yp <- t(matrix(apply(xp, 1, FFfunction, k=i), nrow=ym)) + vks
        f[i, ] <- crossprod(w, yp)
        
        ##a posteriori estimates
        
        #a posteriori state estimate
        m[i + 1, ] <- a[i, ]
        #a posteriori error covariance
        C[[i + 1]] <- R[[i]]
        
        #since there are no measurements available, the importance
        #weights remain unchanged
        
        if (MCparticles) {
          xpt[[i + 1]] <- xp
          wt[i + 1, ] <- w}
      }
      
      else {
        good <- !whereNA
        
        ##time update
        wks <- rmvnorm(n=N, sigma=as.matrix(mod$W))
        #a priori state estimate
        xp <- t(matrix(apply(xp, 1, GGfunction, k=i), nrow=p)) + wks
        a[i, ] <- crossprod(w, xp)
        #covariance of a priori state estimate
        xc <- t(t(xp)-a[i, ])
        R[[i]] <- crossprod(sqrt(w)*xc)
        
        ##measurement update
        vks <- rmvnorm(n=N, sigma=as.matrix(mod$V))
        #predicted measurement
        yp <- t(matrix(apply(xp, 1, FFfunction, k=i), nrow=ym)) + vks
        f[i, ] <- crossprod(w, yp)
        #covariance of predicted measurement
        Qy <- crossprod(sqrt(w)*t(t(yp[, good])-f[i, good]))
        
        ##a posteriori estimates
        
        #importance weights update
        invVgood <- solve(mod$V[good, good])
        xs <- matrix(y[i, good] - t(yp[, good]), nrow=N, byrow=TRUE)
        w <- w*apply(xs, 1, function (x) exp(-.5*(tcrossprod(crossprod((x), invVgood),
                                                             t(x)))))
        w <- w/sum(w)
        
        #a posteriori state estimate
        m[i + 1, ] <- crossprod(w, xp)
        #a posteriori error covariance
        C[[i + 1]] <- crossprod(sqrt(w)*xc)
        
        if (MCparticles) {
          xpt[[i + 1]] <- xp
          wt[i + 1, ] <- w}
        
        ##resampling
        Neff <- 1 / crossprod(w)
        if (Neff < Nthreshold) {
          
          if (resampling=="mult") {
            xp <- as.matrix(resampleM(n=N, wk=w, xk=xp))
            w <- rep(1/N, N)}
          else if (resampling=="sys") {
            xp <- as.matrix(resampleSy(n=N, wk=w, xk=xp))
            w <- rep(1/N, N)}
          else if (resampling=="strat") {
            xp <- as.matrix(resampleSr(n=N, wk=w, xk=xp))
            w <- rep(1/N, N)}
        }
        
        ##roughening
        if (roughening) {
          El <- apply(xp, 2, function (x) diff(range(x)))
          sigma.l <- Grough*El*N^(-1/length(mod$m0))
          if (length(El)==1) Jk <- sigma.l^2 else
            Jk <- diag(sigma.l^2)
          ck <- rmvnorm(n=N, sigma=as.matrix(Jk))
          xp <- xp + ck}
        
        ##compute log-likelihood
        if (logLik) {
          e <- as.matrix(y[i, good] - f[i, good])
          ll <- ll + sum(good)*log(2*pi) + sum(log(eigen(Qy)$values)) + 
            crossprod(e, tcrossprod(solve(Qy, tol=1e-30), t(e)))}
      }
    }
  }
  ans <- list(m = m, C = C, a = a, R = R, f = f)
  
  attributes(ans$f) <- yAttr
  
  if (logLik)
    ans <- c(ans, logLik = 0.5*ll)
  
  if (MCparticles)
    ans <- c(ans, xpt = list(xpt), wt = list(wt))
  
  if (simplify) 
    ans <- c(mod = list(mod1), N = list(N),
             GGfunction = list(GGfunction), FFfunction = list(FFfunction),
             ans)
  else {
    attributes(y) <- yAttr
    ans <- c(y = list(y), mod = list(mod1), N = list(N),
             GGfunction = list(GGfunction), FFfunction = list(FFfunction),
             ans)
  }
  return(ans)
}



resampleM <- function(n, wk, xk) {
  
  #This function implements multinomial resampling
  
  #Arguments:
  #-n is the number of particles
  #-wk is a vector containing the importance weights for the particles
  #-xk is a matrix with the state estimates for the particles
  
  xk <- as.matrix(xk)
  u <- runif(n)
  z <- cumsum(wk)
  i <- sapply(1:n, function (j) sum(z < u[j]))
  xk[i + 1,]
}



resampleSy <- function(n, wk, xk) {
  
  #This function implements systematic resampling
  
  #Last update: 2016/04/13
  
  #Arguments:
  #-n is the number of particles
  #-wk is a vector containing the importance weights for the particles
  #-xk is a matrix with the state estimates for the particles
  
  xk <- as.matrix(xk)
  v <- runif(1)
  u <- sapply(1:n, function (j) ((j-1)+v)/n)
  z <- cumsum(wk)
  i <- sapply(1:n, function (j) sum(z < u[j]))
  xk[i + 1,]
}



resampleSr <- function(n, wk, xk) {
  
  #This function implements stratified resampling
  
  #Last update: 2016/04/13
  
  #Arguments:
  #-n is the number of particles
  #-wk is a vector containing the importance weights for the particles
  #-xk is a matrix with the state estimates for the particles
  
  xk <- as.matrix(xk)
  v <- runif(n)
  u <- sapply(1:n, function (j) ((j-1)+v[j])/n)
  z <- cumsum(wk)
  i <- sapply(1:n, function (j) sum(z < u[j]))
  xk[i + 1,]
}


PFsmooth <- function (filterData) {
  
  #This function implements the backward-simulation particle smoother as
  #described by Godsill, Doucet, and West in their 2004 paper
  #"Monte Carlo Smoothing for Nonlinear Time Series".
  
  #Arguments:
  #-filterdata is an object as returned by PF (with MCparticles=TRUE).
  
  mod <- filterData
  mAttr <- attributes(mod$m)
  mod$m <- as.matrix(mod$m)
  mod$a <- as.matrix(mod$a)
  GGfunction <- mod$GGfunction
  W <- as.matrix(mod$mod$W)
  
  if (all((W)==0)) stop("Backward-simulation is not possible when all entries",
                        "\n  in the W (process noise) matrix are equal to zero")
  invW <- solve(W)
  
  n <- length(mod$R)
  p <- ncol(mod$m)
  
  sT <- mod$xpt[[n + 1]][sample(mod$N, 1, prob=mod$wt[n + 1, ]), ]
  
  s <- rbind(matrix(0, n, p), sT)
  
  
  if (n > 0)
    for (i in n:1) {
      fx <- matrix(apply(mod$xpt[[i]], 1, GGfunction, k=i), nrow=p)
      xs <- t(s[i + 1, ] - fx)
      w <- mod$wt[i, ]*apply(xs, 1, function (x)
        exp(-.5*(tcrossprod(crossprod((x), invW), t(x)))))
      
      s[i, ] <- mod$xpt[[i]][sample(mod$N, 1, prob=w), ]
    }
  
  ans <- list(s = s)
  
  attributes(ans$s) <- mAttr
  return(ans)
}



PFforecast <- function (filterData, nAhead=1) {
  
  #This function predicts future system states and observations for
  #the Particle Filter (PF).
  
  #Arguments:
  #-filterdata is an object as returned by PF (with MCparticles=TRUE).
  #-nAhead is the number of steps ahead for which to predict system states
  # and observations.
  
  mod <- filterData
  
  GGfunction <- mod$GGfunction
  FFfunction <- mod$FFfunction
  N <- mod$N
  
  modFuture <- mod$mod
  nobs <- nrow(as.matrix(mod$a))
  ym <- ncol(mod$f)
  lastObsIndex <- NROW(mod$m)
  modFuture$C0 <- mod$C[[lastObsIndex]]
  modFuture$m0 <- mod$m[lastObsIndex,]
  xp <- mod$xpt[[lastObsIndex]]
  w <- mod$wt[lastObsIndex,]
  
  mod <- modFuture
  p <- length(mod$m0)
  a <- rbind(mod$m0, matrix(0, nAhead, p))
  R <- vector("list", nAhead + 1)
  R[[1]] <- mod$C0
  f <- matrix(0, nAhead, ym)
  Q <- vector("list", nAhead)
  
  for (it in 1:nAhead) {
    
    ##future states
    
    wks <- rmvnorm(n=N, sigma=as.matrix(mod$W))
    #a priori state estimate
    xp <- t(matrix(apply(xp, 1, GGfunction, k=nobs+it), nrow=p)) + wks
    a[it + 1, ] <- crossprod(w, xp)
    #covariance of a priori state estimate
    R[[it + 1]] <- crossprod(sqrt(w)*t(t(xp)-a[it + 1, ]))
    
    ##future observations
    
    vks <- rmvnorm(n=N, sigma=as.matrix(mod$V))
    #predicted measurement
    yp <- t(matrix(apply(xp, 1, FFfunction, k=nobs+it), nrow=ym)) + vks
    f[it, ] <- crossprod(w, yp)
    #covariance of predicted measurement
    Q[[it]] <- crossprod(sqrt(w)*t(t(yp)-f[it, ]))
    
  }
  a <- a[-1, , drop = FALSE]
  R <- R[-1]
  ans <- list(a = a, R = R, f = f, Q = Q)
  
  return(ans)
}