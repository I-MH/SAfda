#  CI.PFF
#
#' Convert the Bootstrapped confidence intervals for the FDF model into 
#' confidence intervals for the PFF model 
#'
#' @param factor_CIs output object from \code{CI.FDF}
#' @param fit.PFF output object from \code{PFF}
#'
#'   A list with: 
#' \itemize{
#'   \item \code{LCIdf} - a data.frame with bootstrapped lower CI for all factors
#'   \item \code{UCIdf} - a data.frame with bootstrapped upper CI for all factors
#' }
#'
#' @examples 
#' library(fda)
#' data(fd.data2011);data(CLRData2011)
#' m<- 2
#' ## Model fit
#' factor_output1 <- FDF(data=fd.data2011,h=24*5,k=m,kmax = 10, 
#'                       kern_type = 'BT')
#' PFF.output<-PFF(factor_output1,dates,CLRData,fd.data)
#' ## CI 
#' set.seed(5643)
#' factor_CIs<-CI.FDF(factor_output1,xarg=CLRData2011$x.fine,n.rep = 10)
#' Sources_CIs<-CI.PFF(factor_CIs=factor_CIs,fit.PFF=PFF.output)
#' 
#' matplot(PFF.output$SourceProfiles, main='Sources',type='l',lty=1)
#' matplot(Sources_CIs$LCIdf[,-1], type='l', lty = 2, add = TRUE)        
#' matplot(Sources_CIs$UCIdf[,-1], type='l', lty = 2, add = TRUE)        

CI.PFF <- function(factor_CIs, fit.PFF) {
  
  #log.x new sequence but can be anything!
  log.x <- fit.PFF$log.x
  
  # num of selected factors 
  kbar <- dim(fit.PFF$SourceProfiles)[2] 
  
  # Transform to Bayes space of factors ------------------------------------------
  hat.fM <- matrix(0, nrow = length(log.x), ncol = kbar)
  Boot.fM <- factor_CIs$Bootstrap.sample
  for (i in 1:kbar) {
    Aux <- matrix(0, nrow = length(log.x),
                  ncol = length(factor_CIs$Bootstrap.sample[1, i, ]))
    for(ib in 1:length(factor_CIs$Bootstrap.sample[1, i, ])) {
      Aux[, ib] <- clr2density(log.x, z_step = fit.PFF$x.step,
                               factor_CIs$Bootstrap.sample[, i, ib])
    }
    Boot.fM[, i, ] <- Aux
  }
  
  #CI ribbons
  LCIdf <- data.frame(x = exp(log.x))
  UCIdf <- data.frame(x = exp(log.x))
  signs <- fit.PFF$signs
  sig2 <- fit.PFF$sig2
  for (j in 1:kbar) {
    Aux <- apply(Boot.fM[, j, ], 2, function(y) {
      prodcf(signs[j] * sig2[j], f = y, s = log.x)$cf })
    df <- data.frame(LCI = fit.PFF$SourceProfiles[, j] - 2 * apply(Aux, 1, sd))
    dfU <- data.frame(UCI = fit.PFF$SourceProfiles[, j] + 2 * apply(Aux, 1, sd))
    LCIdf <- cbind(LCIdf, df)
    UCIdf <- cbind(UCIdf, dfU)
  }
  return(list(LCIdf = LCIdf, UCIdf = UCIdf))
}