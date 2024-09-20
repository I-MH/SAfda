#' PFF
#'
#' After fitting a factor model (output from FDF) in Hilbert space, 
#' it transforms factors and contributions to original density space
#' 
#' @param factor_hat  Output from FDF.R 
#' @param dates Data dates
#' @param CLR_data  Output from PNSDdata.fd
#' @param fd_data   Output from PNSDdata.fd
#' @param flip   tuning parameter to decide if factors should be flipped or not based on 
#' dispersion around main mode. 
#'
#' @return 
#'   A list with: 
#' \itemize{
#'   \item \code{SourceProfiles} - a data.frame with the estimated densities
#'   \item \code{SourceContributions} - a data.frame with dates and time series contributions based on interpretable projection
#'   \item \code{SourceContributions0} - a data.frame with dates and time series contributions based on model fit
#'   \item \code{RowSumSourceContribution} - normalization constant in \code{SourceContributions}
#'   \item \code{log.x, x.step} - values used during density normalization in log range 
#'   \item \code{signs} - (not for users)
#'   \item \code{sig2} - (not for users)
#' }
#'   
#' @export
#'
#' @examples 
#' library(fda)
#' data(fd.data2011);data(CLRData2011)
#' m<- 5
#' factor_output1 <- FDF(data=fd.data2011,h=24*5,k=m,kmax = 10, 
#'                       kern_type = 'BT')
#' PFF.output<-PFF(factor_output1,dates,CLRData,fd.data)
#'                       
#' matplot(PFF.output$SourceProfiles, main='Sources')
#' matplot(PFF.output$SourceContributions[,-1], type='l', 
#'          main=expression(paste('Time series ',beta[t])),
#'          xlab= 'time') 
PFF <- function(factor_hat, dates, CLR_data, fd_data,flip=TRUE){
  
  # log.x new sequence but can be anything!
  log.x <- CLR_data$x.fine
  
  # num of selected factors
  kbar <- dim(factor_hat$hat.F$coefs)[2]   # ****** THIS IS A PROBLEM, NOT DEFINED
  
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
  
  
  if(flip){
    flip.out <- run_flip(SourceProfiles, SourceContributions, log.x)
  }else{
    SourceProfiles<-data.frame(SourceProfiles)
    names(SourceProfiles) <- c(paste0('Factor', 1:kbar))
    flip.out <-list(
      SourceProfiles = SourceProfiles,
      SourceContributions = SourceContributions,
      signs = rep(1,kbar),
      sig2 = NULL
    )
  }
  
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
