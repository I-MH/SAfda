################################################################################
#
#' which_max_v2 An updated `which.max`
#'
#' Computes `which.max`, but remaps boundary values to
#' one-in inset values.
#'
#' @param z Input numeric vector.
#'
#' @return a `which.max` one-in inset from the boundary.
#' @export
#'
#' @examples
#' x <- sample(1:10, 10, TRUE)
#' which_max_v2(x)
#'
which_max_v2 <- function(z) {
  stopifnot(is.numeric(z), length(z) > 2)
  
  nn <- length(z)
  out <- which.max(z)
  if (out == 1)
    out <- 2
  if (out == nn)
    out <- nn - 1
  return(out)
}

################################################################################
#
#' @title clr2density: Invert the clr transformation
#'
#' @param z  Grid of points defining the abscissa
#' @param z_step Step for the grid of the abscissa
#' @param clr   Grid evaluation of the clr transformed density
#'
#' @return  Grid evaluation of the density
#' @export
#'
#' @importFrom fda is.fd
#' @importFrom fda eval.fd
#' 
#' @examples rnorm(1)
clr2density <- function(z, z_step, clr)
{
  if(fda::is.fd(clr)) {
    return(exp(eval.fd(z, clr)) / trapzc(z_step, exp(eval.fd(z, clr))))
  } else {
    return(exp(clr) / trapzc(z_step, exp(clr)))
  }
}

################################################################################
#
#' traczc - Trapezoidal Integration
#'
#' Perform a trapezoidal Riemann sum to the input `y`, with step `step`
#' @param step   Step of the grid.
#' @param y   Vector of points to be integrated.
#'
#' @return  Numerical integration via trapezoidal formula.
#' @export
trapzc <- function(step, y) {
  step * t(y) %*% c(0.5, rep(1, length(y) - 2), 0.5)
}


################################################################################
#' factor_to_sources
#'
#' Function to convert factors back into sources
#' 
#' @param factor_hat   ??
#' @param CLR_data    ??
#' @param fd_data    ??
#'
#' @return list
#' @export
#'
#' @examples rnorm(1)
factor_to_sources <- function(factor_hat, CLR_data, fd_data){

  # log.x new sequence but can be anything!
  log.x <- CLR_data$x.fine
  
  # num of selected factors
  kbar <- m   # ****** THIS IS A PROBLEM, NOT DEFINED

  # Transform to Bayes space of factors
  hat.fM <- matrix(0, nrow = length(log.x), ncol = kbar)
  for (i in 1:kbar) {
    hat.fM[, i] = clr2density(log.x, z_step = CLR_data$x.step,
                              factor_hat$hat.F[i, ])
  }
  
  # Define estimated factors as Source profiles
  SourceProfiles <- hat.fM
  
  # define SourceContributions as the corresponding scores
  SourceContributions <-
    data.frame(dates, factor_hat$hat.beta[, 1:kbar])
  names(SourceContributions) <- c("time", paste0('Factor', 1:kbar))
  
  flip.out <- run.flip(SourceProfiles, SourceContributions, log.x)
  betatilde <- TS.contribution(flip.out$SourceProfiles,
                               flip.out$SourceContributions,
                               fd_data,
                               CLR_data,
                               log.x)

  return(
    list(
      SourceProfiles = flip.out$SourceProfiles,
      SourceContributions = betatilde$SourceContributions3,
      SourceContributions0 = flip.out$SourceContributions,
      RowSumSourceContribution = betatilde$RowSumSourceContribution,
      log.x = log.x,
      signs = flip.out$signs,
      sig2 = flip.out$sig2,
      x.step = CLR_data$x.step
    )
  )
}

################################################################################
#' run_flip
#'
#' @param SourceProfiles fill in
#' @param SourceContributions fill in
#' @param log.x fill in
#'
#' @return list
#' @export
run_flip <- function(SourceProfiles,
                     SourceContributions,
                     log.x) {

  kbar <- ncol(SourceProfiles) # num of selected factors

  # Flip only considering profiles
  # The criteria takes the current SP and looks for the flip that shows a max mode
  # preference to concave shape like than convex?
  flip <- function(SP, der.tol = 0.1) {
    op1 <- SP
    op2 <- prodcf(-1, f = SP, s = log.x)$cf
    pos.m1 <- which.max.v2(op1)
    pos.m2 <- which.max.v2(op2)

    der1 <-
      abs((op1[pos.m1 + 1] - op1[pos.m1 - 1]) / (log.x[pos.m1 + 1] - log.x[pos.m1 -
                                                                             1]))
    der2 <-
      abs((op2[pos.m2 + 1] - op2[pos.m2 - 1]) / (log.x[pos.m2 + 1] - log.x[pos.m2 -
                                                                             1]))

    if (der1 < der.tol & der2 < der.tol) {
      out.sign <- c(1, -1)[which.max(c(op1[pos.m1], op2[pos.m2]))]
    } else {
      out.sign <- c(1, -1)[which.min(c(der1, der2))]
    }
    return(out.sign)
  }
  signs0 <- apply(SourceProfiles, 2, FUN = flip)

  #Correct Flip on positive betas betas
  #require as negative always beta means a flipped true source
  SourceContributions2 <- SourceContributions
  for (j in 1:kbar) {
    SourceContributions2[, j + 1] <-
      SourceContributions[, j + 1] * signs0[j]
  }

  change.sign <- which(colMeans(SourceContributions2[, 2:(kbar + 1)]) < 0)
  signs1 <- signs0
  signs1[change.sign] <- (-1) * signs0[change.sign]
  signs <- signs1

  #Betas with (mean one)
  #require to minimize dispersion around mode (if any!)
  SourceContributions2 <- SourceContributions
  for (j in 1:kbar) {
    SourceContributions2[, j + 1] <- SourceContributions[, j + 1] * signs[j]
  }
  sig2 <- rep(1, kbar)
  means.aux <- colMeans(SourceContributions2[, 2:(kbar + 1)])
  sig2[means.aux > 0.5] <-
    colMeans(SourceContributions2[, 2:(kbar + 1)])[means.aux > 0.5]
  Factdf <- data.frame(x = exp(log.x))
  for (j in 1:kbar) {
    var <- paste0('Factor', j)
    df <-
      data.frame(Factor = prodcf(signs[j] * sig2[j], f = SourceProfiles[, j], s =
                                   log.x)$cf)
    Factdf <- cbind(Factdf, df)
  }
  nameF <- paste0('Factor ', 1:kbar)
  names(Factdf) <- c('x', paste0('Factor', 1:kbar))
  SourceProfiles2 <- Factdf[, -1]
  # fix signs for source contribution accordingly
  SourceContributions2 <- SourceContributions
  for (j in 1:kbar) {
    SourceContributions2[, j + 1] <-
      SourceContributions[, j + 1] * signs[j] * (1 / sig2[j])
  }
  return(
    list(
      SourceProfiles = SourceProfiles2,
      SourceContributions = SourceContributions2,
      signs = signs,
      sig2 = sig2
    )
  )
}

################################################################################
#' TS.contribution
#' 
#' Blah blah blah
#'
#' @param SourceProfiles2 fill in
#' @param SourceContributions2 fill in
#' @param fd.data fill in
#' @param CLR_data fill in
#' @param log.x fill in
#'
#' @return fill in
#' @export
#'
#' @examples rnorm(1)
TS.contribution <-
  function(SourceProfiles2,
           SourceContributions2,
           fd.data,
           CLR_data,
           log.x) {
    
    # num of factors
    nprofiles <- dim(SourceProfiles2)[2] 
    
    # Data as Functional object
    # need fd.data
    Clr.fdfit <- fda::eval.fd(evalarg = CLR_data$x.fine, fdobj = fd.data)
    fdfit <-
      apply(Clr.fdfit,
            2,
            clr2density,
            z = CLR_data$x.fine,
            z_step = CLR_data$x.step)

    varNames <- names(SourceContributions2)[-1]
    SourceContributions3 <- data.frame(time = SourceContributions2$time)
    delta <- diff(log.x)[1]
    
    for(k in 1:nprofiles) {
      f1Xi <- apply(fdfit, 2, function(x)
          x * SourceProfiles2[, k])   # X_i * f_1
      SourceContributions3[varNames[k]] <- delta * apply(f1Xi, 2, sum)
    }
    
    # normalize contributions such that sums to one over the diff sources
    RowSumSourceContribution3 <-  rowSums(SourceContributions3[, -1])
    SourceContributions3[, -1] = t(apply(SourceContributions3[, -1], 1, function(x)
      x / sum(x)))
    return(
      list(
        SourceContributions3 = SourceContributions3,
        RowSumSourceContribution = RowSumSourceContribution3
      )
    )
  }

################################################################################
#  FactorCIToSourcesCI
#
#' Convert Factor CIs to Source CIs
#'
#' @param factor_CIs fill in
#' @param Res.FTS fill in
#'
#' @return list
#' @export
#'
#' @examples rnorm(1)
FactorCIToSourcesCI <- function(factor_CIs, Res.FTS) {

  #log.x new sequence but can be anything!
  log.x <- Res.FTS$log.x
  
  # num of selected factors 
  kbar <- dim(Res.FTS$SourceProfiles)[2] 
  
  # Transform to Bayes space of factors ------------------------------------------
  hat.fM <- matrix(0, nrow = length(log.x), ncol = kbar)
  Boot.fM <- factor_CIs$Bootstrap.sample
  for (i in 1:kbar) {
    Aux <- matrix(0, nrow = length(log.x),
                  ncol = length(factor_CIs$Bootstrap.sample[1, i, ]))
    for(ib in 1:length(factor_CIs$Bootstrap.sample[1, i, ])) {
      Aux[, ib] <- clr2density(log.x, z_step = Res.FTS$x.step,
                               factor_CIs$Bootstrap.sample[, i, ib])
    }
    Boot.fM[, i, ] <- Aux
  }
  
  #CI ribbons
  LCIdf <- data.frame(x = exp(log.x))
  UCIdf <- data.frame(x = exp(log.x))
  signs <- Res.FTS$signs
  sig2 <- Res.FTS$sig2
  for (j in 1:kbar) {
    Aux <- apply(Boot.fM[, j, ], 2, function(y) {
        prodcf(signs[j] * sig2[j], f = y, s = log.x)$cf })
    df <- data.frame(LCI = Res.FTS$SourceProfiles[, j] - 2 * apply(Aux, 1, sd))
    dfU <- data.frame(UCI = Res.FTS$SourceProfiles[, j] + 2 * apply(Aux, 1, sd))
    LCIdf <- cbind(LCIdf, df)
    UCIdf <- cbind(UCIdf, dfU)
  }
  return(list(LCIdf = LCIdf, UCIdf = UCIdf))
}

################################################################################
#' data.fd
#'
#' Form a `fd` object from input data.
#' 
#' @param argvals fill in
#' @param data fill in
#' @param basis fill in
#' @param nbasis fill in
#' @param rangeval fill in
#' @param lambda fill in
#'
#' @return fill in
#' @export
#'
#' @examples rnorm(1)
data.fd <- function(argvals = NULL,
                    data,
                    basis,
                    nbasis = 25,
                    rangeval = c(0, 1),
                    lambda = NULL) {
  # returns data with class fd
  if (is.numeric(data)) {
    N <- dim(data)[2]
    m <- dim(data)[1]
    datafd0 <- SmoothData(
      argvals = argvals,
      data = data,
      type_basis = basis,
      nbasis = nbasis,
      rangeval = rangeval,
      lambda = lambda
    )$fd
  } else{
    if (sum(class(data) %in% "fd") != 0) {
      rangeval = data$basis$rangeval
      N <- dim(data$coefs)[2]
      datafd0 <- data
      m <- 100
      nbasis <- data$basis$nbasis
    } else{
      stop("data should be a matrix or fd object")
    }
  }
  return(list(
    datafd = datafd0,
    rangeval = rangeval,
    N = N,
    m = m,
    nbasis = nbasis
  ))
}

################################################################################
#' diff_fd
#'
#' 1-step difference of a data set in class `fd` format, returned as `fd`.
#'
#' @param datafd Input data in FD format.
#'
#' @importFrom fda fd
#' @return Differentiated (1-step differencing operator) input, in `fd` format.
#' @export
#'
#' @examples rnorm(1)
diff_fd <- function(datafd) {

  # it computes X_n - X_{n-1}
  # datafd= data of class fd from data.fd
  newCoeff <- t(apply(datafd$coefs, 1, diff))
  new.fd <- fda::fd(newCoeff, datafd$basis)
  new.N <- dim(newCoeff)[2]
  datafd$coefs <- newCoeff
  return(datafd)
}

# Bayes operations
# f g assume to be discrete

#' sumfg
#'
#' Utility function to compute the midpoint-rule Riemann sum of `f * g`
#' with possibly irregular bases (from the first-difference of `s`) and 
#' then return the product of `f * g`, normalized by that sum.
#'
#' @param f Input 1
#' @param g Input 2
#' @param s Observation times; first difference is the base for the sum.
#'
#' @return A list of two objects: the `sum` and the `norm` (constant)
#' @export
#'
#' @examples ss<-seq(-4,5,length.out=100)
#' f<-dnorm(ss,mean=0,sd=1)
#' g<-dnorm(ss,mean=1,sd=sqrt(2))
#' fpg<-sumfg(f,g,ss)$sum
#' plot(ss,f,type="l",lwd=2,col="orange",ylim=c(0,max(f,g,fpg)),ylab="",main="f plus g")
#' lines(ss,g,lwd=2,col="blue")
#' lines(ss,fpg,lwd=2,col="green")
#' legend('topright',legen=c("f","g"),col=c("orange","blue"),lty=1)

sumfg <- function(f,g,s){
  stopifnot(is.numeric(f), is.numeric(g), is.numeric(s)) 
  # Assuming s is sorted - check
  stopifnot(all(diff(s) > 0))
  
  aux.fg <- f * g
  n <- length(s)
  # int <- (s[2:n]-s[1:(n-1)]) %*% (aux.fg[2:n]+aux.fg[1:(n-1)])/2)
  int <- diff(s) %*% ((aux.fg[2:n] + aux.fg[1:(n-1)]) / 2)
  return(
    list(
      sum = aux.fg / c(int), 
      norm = int
    )
  )
}

#' sumF
#'
#' Act across the columns of `Fmat`, and compute the `sumfg` on
#' the columns pairwise, starting from the first two. Continue
#' doing so until all columns are multiplied. Normalizes each product.
#'
#' @param Fmat 
#' @param s Observation times; first difference is the base for the sum.
#'
#' @return The normalized vector of products.
#' @export
#'
#' @examples rnorm(1)
sumF <- function(Fmat, s){
  nf <- ncol(Fmat)
  if(nf == 2) {
    result <- sumfg(Fmat[, 1], Fmat[, 2], s)$sum 
  } else { 
    tempf <- Fmat[, -nf]
    # recurse on columns 1:(nf-1)
    tempf2 <- sumF(tempf, s)
    result <- sumfg(Fmat[, nf], tempf2, s)$sum
  }
  return(result)
}

#' prodcf
#'
#' Utility function to compute the midpoint-rule Riemann sum of `f^c`
#' with possibly irregular bases (from the first-difference of `s`) and 
#' then return the product of `f^c`, normalized by that sum.
#'#'
#' @param c Power
#' @param f Base
#' @param s Observation times; first difference is the base for the sum.
#'
#' @return A list of two objects: the `sum` and the `norm` (constant)
#' @export
#'
#' @examples ss<-seq(-4,5,length.out=100) 
#' f<-dnorm(ss,mean=0,sd=1) alp<-2
#' alp.times.f<-prodcf(alp,f=f,s=ss)$cf
#' plot(ss,f,type="l",lwd=2,col="orange",ylim=c(0,max(f,alp.times.f)))
#' lines(ss,alp.times.f,lwd=2,col="green")
prodcf <- function(c, f, s) {
  # Assuming s is sorted
  stopifnot(all(diff(s) > 0))
  stopifnot(is.numeric(c), is.numeric(f), is.numeric(s))
  
  aux.cf <- f^c 
  n <- length(s)
  # int <- as.numeric((s[2:n]-s[1:(n-1)])%*%(aux.cf[2:n]+aux.cf[1:(n-1)])/2)
  int <- diff(s) %*% (aux.cf[2:n] + aux.cf[1:(n - 1)]) / 2
  return(
    list(
      cf = aux.cf / c(int),
      norm = int
    )
  )
}

