# ncodages.dat from steve
# get into right dimensions, feed to plotWt.R
# plots alpha and rho over time, writes a file
# usually not well estimated for the last few cohorts, go back 5-10 years
# and copy stable estimates forward
# fit mseR's age structured model to cod data: populate 
# assessCA.dat with data from steve. Only need fishery catch series,
# 1 index series, there'll be some matching
# Fishery age comp, survey age comp
# Can growth pars from plotWt, as long as we have the c1/c2 pars
# c2 is length vs weight relationship (power)
# c1 is kg to cm
# estimate til later (til we get it right)
# assessCA.r will grab rep file and plot, so run ADMB separately.

c1 <- 0.0000741
c2 <- 3.06
makeSizeAge <- function( fileName="weightAgeSurvey.txt", yrs=c(1959,2014) )
{
  tmp <- read.table(file=fileName,sep="\t")
  Wta <- as.matrix( tmp[,-1] )

  ages   <- 1:ncol(Wta)
  maxWt  <- max(Wta,na.rm=TRUE)
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )

  dim <- yrs [ 2 ] - yrs [ 1 ] + 1

  matYr <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matWt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matLt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  for( cohRow in 1:dim )
  {
    maxCol <- min( cohRow + maxAge - 1, ncol(matWt) )
    a <- 0
    for( j in cohRow:maxCol )
    {
      a <- a + 1
      matYr[ cohRow, j ] <- (yrs[1] + cohRow -1) + a - 1
      matWt[ cohRow, j ] <- Wta[ cohRow, a ]
      logL               <- (log(Wta[cohRow,a])-log(c1))/c2
      matLt[ cohRow, j ] <- exp(logL)
    }
  }
  result <- list()
  result$Yr <- matYr
  result$Wt <- matWt
  result$Lt <- matLt
  return( result )

}  
plotMat <- function( matYr, mat, yrs=c(1959,2014) )
{
  ages   <- 1:ncol(mat)
  maxX   <- max(mat,na.rm=TRUE)
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )
   xLim  <- yrs
   yLim  <- c( 0,maxX )
   pCol  <- rep( c("white","white","white","white","black"), 8 )

  # Weight at age.
  plot( xLim,yLim,type="n",axes=FALSE,xlab="",ylab="" ) 
  
  for( t in 1:nrow(mat) )
  {
    lines( matYr[t,], mat[t,], lty=1, lwd=1 )
    points( matYr[t,], mat[t,], bg=pCol[t], col="black", cex=.8, pch=21 )
  }

  axis( side=1, cex.axis=1 )
  axis( side=2, cex.axis=1, las=1 )
  box()
} 
fitWalford <- function( mat, yrs=c(1970,2011) )
{
  dim <- yrs [ 2 ] - yrs [ 1 ] + 1
  coeffs <- matrix(NA,nrow=dim,ncol=3)
  len    <- ncol(mat)
  yrs    <- yrs[1]:yrs[2]
  for( t in 1:nrow(mat) )
  {
    la1 <- mat[t,1:(len-1)]
    la2 <- mat[t,2:len]
    if( t < nrow(mat) )
    {
      tmp <- lm(la2~la1)
      coeffs[t,1]   <- yrs[t]
      coeffs[t,2:3] <- tmp$coef[1:2]
    }
  }
  return( coeffs )
}



# par( oma=c(3.5,4,1,1), mar=c(2,2,1,1), mfrow=c(2,1 ) )
sizeListS <- makeSizeAge( fileName="weightAgeSurvey.txt", yrs=c(1959,2014) )
# plotMat(matYr=sizeListS$Yr, mat=sizeListS$Wt)
# legend(x=1970,y=6,legend="Survey",bty="n")
# # sizeList <- makeSizeAge( fileName="weightAgeFishery.txt", yrs=c(1959,2014) )
# # plotMat(matYr=sizeList$Yr, mat=sizeList$Wt)
# # legend(x=1970,y=6,legend="Fishery",bty="n")
# # mtext( side=1, line=1, cex=1, outer=TRUE, "Cohort/Year" )
# # mtext( side=2, line=1, cex=1, outer=TRUE, "Weight (kg)" )

walPars <- fitWalford(mat=sizeListS$Lt, yrs=c(1959,2014))

# plot ( 1959:2014, y = walPars[,2] )

