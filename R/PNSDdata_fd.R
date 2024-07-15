PNSDdata_fd<-function(Dens,x,n.basis=20,k = 4,
                      der = 2,
                      alpha = 0.999,
                      ch = 1 ){
  #Dens: matrix per row has each PNSD. Dens[i,j] gives the PNSD at time i and particle size j
  #x: vector of particle sizes  
  #         k = order of spline, degree = k-1
  #         der = derivation.
  #         alfa = smoothing parameter
  #         ch ={1,2};  ch=1: functional with (1-alfa) and alfa
  #                     ch=2: funkcional with alfa   
  
  #Define grid for x in log scale------------------------------------------------
  # If data is integer then lowbound.x has to be changed! to x-epsilon
  lowbound.x<-floor(x)[1]
  log.x=c(log(lowbound.x),log(x))
  # x.fine must be of length 1000 to match SmoothingSpline.R
  x.fine = seq(min(log.x), max(log.x), length=1000) 
  x.step = diff(x.fine[1:2])
  ### widths + volumes of intervals in histogram
  width=as.matrix(diff(log.x)) 
  centers=NULL
  for (i in 1:length(log.x)-1){
    centers[i]=log.x[i]+(log.x[i+1]-log.x[i])/2
  }
  centers=as.matrix(centers) 
  
  #Normalize PNSD (Dens matrix)--------------------------------------------------
  norm.hist=matrix(nrow=nrow(Dens),ncol=ncol(Dens))
  c<-numeric(nrow(Dens))
  for(i in 1:nrow(Dens)){
    c[i]<-sum(Dens[i,])
    norm.hist[i,]=as.numeric(Dens[i,]/c[i]) 
  }
  densities=matrix(nrow=nrow(Dens),ncol=ncol(Dens))
  for (j in 1:ncol(densities)){
    densities[,j]=norm.hist[,j]/width[j]
  }
  
  #Clr transformation using robCompositions package -----------------------------
  Dens.clr = robCompositions::cenLR((densities))$x.clr
  
  #B-splines parameters-----------------------------------------------------------
  #n.basis=20
  knots=seq(min(x.fine),max(x.fine),length=n.basis) 
  w = rep(1,ncol(Dens.clr)) 
  k = 4
  der = 2
  alpha = 0.999
  ch = 1     
  t =c(t(centers))
  # inputs:  knots= knots of the spline
  #         t = point of approximation,
  #         f = values at t,
  #         w = coefficients for weights,
  #         k = order of spline, degree = k-1
  #         der = derivation.
  #         alfa = smoothing parameter
  #         ch ={1,2};  ch=1: functional with (1-alfa) and alfa
  #                     ch=2: funkcional with alfa   
  
  #generating z-coeficients for B-spline basis------------------------------------
  J=nrow(Dens.clr)
  z_coef=matrix(nrow=length(knots)+1,ncol=nrow(Dens.clr))
  for (i in 1:J){
    # SmoothingSpline0 defined in SmoothingSpline.R
    AuxSplinefit<-smoothing_spline0(knots=knots, tp=t, f=as.numeric(Dens.clr[i,]), w=w, 
                                    k=k, der=der, alpha=alpha, ch=ch)
    z_coef[,i]=AuxSplinefit$z
    J[i]=AuxSplinefit$J
  }
  z_coef=as.matrix(z_coef)
  ###### ZB-spline basis  ######
  #ZsplineBasis defined in NulBaze.R
  AuxZB<-z_spline_basis(knots = knots,k)
  Z = AuxZB$C0
  b_coef = t(AuxZB$D)%*%AuxZB$K%*%z_coef
  
  #Output arrangement -----------------------------------------------------------
  #Transformed data
  CLRData<-vector("list")
  CLRData$x.fine<-x.fine
  CLRData$x.step<-x.step
  CLRData$centers<-centers
  CLRData$densities<-densities
  CLRData$Dens.clr<-Dens.clr
  #Basis info
  Bsplinepar<-vector("list")
  Bsplinepar$knots<-knots
  Bsplinepar$w <- w
  Bsplinepar$k <- k
  Bsplinepar$der <- der
  Bsplinepar$alfa <- alpha
  Bsplinepar$ch <- ch     
  Bsplinepar$t <- t
  #fd.data
  # B-spline basis -------------------------------------------------------
  B = create.bspline.basis(range(knots), nbasis = dim(z_coef)[1]+1,
                           norder=k,breaks=knots)
  # Data as Functional object-----------------------------------------------------
  fd.data = fd(b_coef,B) 
  return(list(z_coef,Z,b_coef,Bsplinepar,CLRData,fd.data))
}