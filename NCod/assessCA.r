lninvgamm<-function(x,alpha=1,beta=1)
{
  logpdf<-alpha*log(beta)-(alpha+1)*log(x)-log(gamma(alpha))-beta/x
  return( logpdf )
}
  
lisread <- function( fname,quiet=FALSE )
{
  # lisread: Function to read a list of data objects from a file.
  # The initial characters "##" denote a comment line (ignored).
  # The initial characters "# " denote a variable name.
  # All other lines must contain scalars or vectors of numbers.
  # Furthermore, all rows in a matrix must contain the same number of
  # columns. Row and column vectors are not converted to matrices.
  #
  # fname  : File name.
  # quiet  : If true, shut up about reporting progress.
  # result : List object with components named in the file.

  # Original functions courtesy of Jon Schnute.
  # Modifications by A.R. Kronlund.

  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x,white=T )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )
    xc
  }

  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables
  nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1     # number of rows in var i
  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1
    i2 <- i1 + nrow[i] - 1                       # range of lines for var i
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}

calcSteepnessNLP <- function( logit.ySteepness=0.5, pmSteepness=0.6, pvSteepness=0.005 ){

  ySteepness <- exp(logit.ySteepness)/(1.+exp(logit.ySteepness))
  rSteepness <- (ySteepness+0.25)/1.25

  # Beta prior on steepness
  muB <- 1.25*pmSteepness-0.25
  tauB <- muB*(1.-muB)/(1.5625*pvSteepness)-1.
  aB <- tauB*muB
  bB <- tauB*(1.-muB)
  steepnessNLP = (-1.)*( (aB-1.)*log(ySteepness) + (bB-1.)*log(1.-ySteepness) )
  
  return( steepnessNLP )

}

plotRtFt <- function( fitObj )
{


  omega <- fitObj$omega
  Ft    <- fitObj$Ft

  lenOmega <- length(omega)
  lenFt    <- length(Ft)
  
  par(mfrow=c(2,1), mar=c(4,3,1,1), mar=c(4,4,1,1))
  
  plot( 1960+1:lenOmega, omega, type="b", pch=19, ylim=c(-1,1),
        ylab="log-recruitment deviation",
        xlab="Year")
  legend(x=0,y=1,legend="A",bty="n")
        
  plot( 1960+1:lenFt, Ft, type="b", pch=19, ylim=c(0,0.8),
        xlab="Year", ylab="Fishing mortality rate")  
  legend(x=0,y=.8,legend="B",bty="n")  
}


result <- lisread("assessca.rep", quiet=TRUE)


plotAll <- function( fitObj )
{
  nT <- result$nT
  nAges <- result$nAges
  par( oma=c(2,2,1,1), mar=c(4,4,1,1), mfcol=c(3,2) )

  result$ItScaled[ result$ItScaled < 0 ] <- NA

  yLim <- c(0,max( max(na.omit(result$ItScaled)), max(result$exploitBt) ) )
  plot( c(1:nT), result$exploitBtS, type="l", lwd=2, ylim=yLim )
    points( c(1:nT), result$ItScaled )
    lines( c(1:nT), result$ItScaled )
    lines( c(1:nT), result$spawnBt, col="blue" )

  yMax <- max(result$Rt)
  plot( c(1:nT), result$Rt[1:nT], type="b", ylim=c(0,yMax) )

  B <- seq(0,result$B0,length=100)
  rec <- result$rec.a*B/( 1.0 + result$rec.b*B )
  plot( B,rec, type="l", ylim=c(0,max(result$Rt) ) )
    points( result$spawnBt[1:nT], result$Rt[1:nT] )

  plot( x=NULL, y=NULL, xlim=c(1,nT), ylim=c(1,nAges),
        xlab="Year",ylab="Age", main="obsPropAge" )
  for( t in 1:nT )
  {
    points( rep(t,nAges),c(1:nAges), cex=10.*result$obsPropAge[t,] )	

  }

  plot( x=NULL, y=NULL, xlim=c(1,nT), ylim=c(1,nAges),
        xlab="Year",ylab="Age", main="predPropAge" )
  for( t in 1:nT )
  {
    points( rep(t,nAges),c(1:nAges), cex=10.*result$predPropAge[t,] )	

  }


  plot( x=NULL, y=NULL, xlim=c(1,nT), ylim=c(1,nAges),
        xlab="Year",ylab="Age", main="residPropAge" )
  resCol <- vector("numeric")
  for( t in 1:nT )
  {
  	resid <- 5.*(result$obsPropAge[t,]-result$predPropAge[t,])
  	resCol[resid > 0] <- "blue"
  	resCol[resid <= 0] <- "red"
  	points( rep(t,nAges),c(1:nAges), cex=resid, col=resCol )	
  }
  
}



