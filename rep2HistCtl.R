# --------------------------------------------------------------------
# rep2HistCtl.R
#
# A script to transform a given ADMB rep file for a stock into
# a simhistory file and a base simCtlFile.
#
# --------------------------------------------------------------------

# Clean house
# rm ( list = ls () )

# Load packages
source ( "msertools.R" )

# Function to create simHistory file from a rep file
# containing recruitment deviations, mortality and 
# numbers at age. Also requires weight at age data from
# another dat file
              ###### NOTE ######
# Requires estimates of c1 and c2 for weight-length relationship
# W = c1 * L^{c2}, where c1 is converting cm to kg, not mm to g.

##### Enter values below before running script #####

# weight-at-length parameters - taken from fishbase
c1 <- 0.00000741                # unit conversion (cm to kg)
c2 <- 3.06                      # "Volume" conversion

# Data files in working directory
repfile <- "ncod.rep"           # rep file from operating model
agefile <- "ncodAges.dat"       # age dat file, containing weight at age data
idxfile <- "ncodIndex_sim.dat"  # dat file containing indices of abundance
katfile <- "ncodCatch.dat"      # catch dat file

# Set time extent of rep file and first year of history file
tExt <- c(1959,2014)

# Base simCtl file in working directory
ctlFile <- "simCtlFile.txt"


##### Functions for history file ######

# This function makes size, age and cohort matrices, used
# in fitWalford to create Walford intercepts and slopes for
# a simHistory file.
makeSizeAge <- function( wtAgeMat = t(wtAge), yrs=c(1959,2014) )
{
  # remove first column of weight at age matrix (first age class)
  Wta <- as.matrix( wtAgeMat[,-1] )

  # Recover ages, max age, max weight
  ages   <- 1:ncol(Wta)
  maxWt  <- max(Wta,na.rm=TRUE)
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )

  # Compute the dimension of Year, Weight and Length matrices
  dim <- yrs [ 2 ] - yrs [ 1 ] + 1

  # Initialise matrices
  matYr <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matWt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )
  matLt <- matrix(NA, nrow=dim, ncol=dim, byrow=T )

  # Loop to populate matrices
  for( cohRow in 1:dim )
  {
    maxCol <- min( cohRow + maxAge - 1, ncol(matWt) )
    a <- 0
    for( j in cohRow:maxCol )
    {
      # Fill cohort, weight and length matrices using 
      # weight at age matrix and weight at length pars.
      a <- a + 1
      matYr[ cohRow, j ] <- (yrs[1] + cohRow -1) + a - 1
      matWt[ cohRow, j ] <- Wta[ cohRow, a ]
      logL               <- (log(Wta[cohRow,a])-log(c1))/c2
      matLt[ cohRow, j ] <- exp(logL)
    }
  }
  # Return matrices.
  result <- list()
  result$Yr <- matYr
  result$Wt <- matWt
  result$Lt <- matLt
  return( result )
}  

# PlotMat will plot the entries in either a weight or length
# at age matrix, built by makeSizeAge()
plotMat <- function( matYr, mat, yrs=c(1959,2014) )
{
  # Recover ages, ranges and dimensions of matrices.
  ages   <- 1:ncol(mat)
  maxX   <- max(mat,na.rm=TRUE)
  maxAge <- max(ages)
  xRange <- c( 0,maxAge )
   xLim  <- yrs
   yLim  <- c( 0,maxX )
   pCol  <- rep( c("white","white","white","white","black"), 8 )

  # Make an empty plot.
  plot( xLim,yLim,type="n",axes=FALSE,xlab="",ylab="" ) 
  
  # Now plot the entries in the given matrix
  for( t in 1:nrow(mat) )
  {
    lines( matYr[t,], mat[t,], lty=1, lwd=1 )
    points( matYr[t,], mat[t,], bg=pCol[t], col="black", cex=.8, pch=21 )
  }

  # Place axes
  axis( side=1, cex.axis=1 )
  axis( side=2, cex.axis=1, las=1 )
  box()
} 

# fitWalford()
# A function that computes the Walford intercepts and slopes for each year
# by fitting a linear model between one year's length at age and
# the next for every age class in a given year. Due to the shrinking length of 
# time series in the final 20 years, estimates tend to destabilise
# towards the end of the time series.
# inputs:     mat = matrix of length at age produced by makeSizeAge()
#             yrs = numeric of min and max years of time series.
# outputs:    coeffs = matrix with walford int and slope for all but 2 years
# usage:      Walford parameters are used in mseR sim history files for 
#             individual growth
fitWalford <- function( mat, yrs=c(1970,2011) )
{
  # Compute dimension of output
  dim <- yrs [ 2 ] - yrs [ 1 ] + 1

  # Initialise coefficients matrix for output
  coeffs <- matrix(NA,nrow=dim,ncol=3)
  
  len    <- ncol(mat)      # Get the length of the series
  yrs    <- yrs[1]:yrs[2]  # create a vector of years
  
  # Loop to fill walford coefficients
  for( t in 1:nrow(mat) )
  {
    # Get time series of length at age, offset by 1
    la1 <- mat[t,1:(len-1)]
    la2 <- mat[t,2:len]
    if( t < nrow(mat) )
    {
      # Fit a linear model
      tmp <- lm(la2~la1)
      # reocver coefficients and populate output matrix
      coeffs[t,1]   <- yrs[t]
      coeffs[t,2:3] <- tmp$coef[1:2]
    }
  }
  return( coeffs )
}

# histCreate()
# A function to write a simhistory file built from a given
# ADMB rep file and age/weight structure dat file. 
# inputs:     repFile = stock's ADMB report file;
#             ageFile = stock's ADMB input age structure
#             yrs = vector of min/max years; 
#             fyr = first year of output history file 
#                   (if NULL defaults to min year)
#             outFile = character vector of output file name
# output:     None
# usage:      Used to write a mseR simHistory file.
histCreate <- function (  rep = lisread ( repfile ), 
                          ages = lisread ( agefile ),
                          yrs = c(1959,2014), fyr = NULL,
                          outFile = "simHistNCod83.csv" )
{
  # Compute length of simHistory file
  nYrs <- yrs [ 2 ] - yrs [ 1 ] + 1

  # Set start year, defaults to the min year unless fyr provided
  start <- 1
  if ( !is.null ( fyr ) ) 
  {
    if ( fyr < yrs [ 1 ] ) cat ("Provided fyr is too early, try again. \n")
    start <- fyr - yrs [ 1 ] + 1
  }

  # recover time series length
  nF <- length ( rep $ Ftg )
  if ( nF != nYrs ) 
  {
    cat ( "Dimension mismatch, check year bounds and try again.\n" )
    return()
  }

  # Now start recovering the "straightforward" values
  Ft     <- rep $ Ftg
  Mt     <- rep $ Mt[2,]
  omegat <- rep $ recDevs

  # Pad Omega (shorter than Ft)
  omegat <- c(0, omegat)

  # Now do growth parameters from fitWalford()
  wtAge <- t ( ages $ wtAge )

  # replace the first few years with the averages for
  # the remaining years
  # Set max year for replacement
  repYear <- 25 
  # Compute averages
  aveWt <- apply ( X = wtAge [ , (repYear+1):ncol(wtAge) ], 
                   FUN = mean, MARGIN = 1 )
  # Replace early years with later year averages.
  wtAge [ , 1:repYear] <- aveWt

  # Create list of matrices of size, age, cohorts.
  sizeListS <- makeSizeAge( wtAgeMat = t(wtAge), yrs=yrs )
  walPars <- fitWalford(mat=sizeListS$Lt, yrs=yrs)

  # Extract walford intercept and slope
  alphat <- walPars [ 1:(nYrs-10), 2 ]
  rhot   <- walPars [ 1:(nYrs-10), 3 ]

  # Replace wild final years of alphat and rhot
  alphat[(nYrs-10 + 1):nYrs]  <- alphat[(nYrs-20 + 1):(nYrs-10)]
  rhot[(nYrs-10 + 1):nYrs]    <- rhot[(nYrs-20 + 1):(nYrs-10)]

  # Now calculate initN - multipliers to convert initial
  # equilibrium abundance to given numbers at age
  R0      <- rep $ R0                 # Recover equilibrium recruitment
  M       <- ( rep $ M [ 1 ] )        # Recover mortality
  Nt      <- rep $ Nta[start,]        # Get numbers at age for start year
  Neq     <- R0 * exp ( - M * 0:19 )  # Compute equilibrium age structure
  initN   <- Nt / Neq                 # Compute initN mults

  # Pad out initN with NAs
  nPad  <- max ( 0, nYrs - start + 1 - 20 )
  pad   <- rep ( NA, nPad )
  initN <- c( initN, pad )

  # Create a data.frame to hold the info
  simHist <- data.frame ( omegat = omegat [ start:nYrs],
                          Ft = Ft [ start:nYrs],
                          Mt = Mt [ start:nYrs],
                          alphat = alphat [ start:nYrs],
                          rhot = rhot [ start:nYrs],
                          initN = initN )

  # Write to file.
  cat ( "Writing to ", outFile, "\n", sep = "" )
  write.csv ( file = outFile, x = simHist )              
  return(simHist)
}

##### Now to create a simCtlFile #####

## things to update in the control file (opMod):
# 1. tMP = length history  + 1
# 2. nT = tMP + 29
# 3. history file
# 4. growth pars:
#   - vonK, Linf, L1
#   - c1, c2
# 5. Selectivity and Maturity
# 6. nAges
# 7. B0
# 8. Steepness

# MP changes:
# 1. t1method(=tMP), t2method(=nT)
# 2. t2Ages, t2AgesS, t2Survey
# 3. FtEndYr
# 4. upperBaseEndYr, lowBaseEndYr

# The following functions are used to create the control list, which
# is then saved out to the filename specified as simCtlName at the top
# of this script.


# saveSimPars   (saves the Sim GUI parameters to a file for runMSE loop).
# Purpose:      Saves an ASCII file with simulation controle parameters.
# Parameters:   Output file name for ASCII file containing GUI parameters.
# Returns:      NULL
# Side effects: Writes an ASCII control file configured for read into a list.
# Source:       A.R. Kronlund
saveSimPars <- function( ctlPars, parFile="inputParameters.txt", overWrite=FALSE )
{
  if ( !overWrite )
  {
    tmpFile <- selectFile( initialfile=parFile,
                  filetype=list( c(".txt","txt files")), mode="save" )

    if ( is.null(tmpFile) )
      return( invisible() )
    else
      parFile <- tmpFile
  }
  else
  {
    # Force an overwrite of parFile for runMSE.
    if ( .WHINE > 1 )
      cat( "\nMSG (saveSimPars) Overwriting ",parFile,"\n" )
  }

  cat( file=parFile, paste( "# mseRguiSim GUI parameters written",date(),"\n" ) )
  cat( file=parFile, "parameter value\n", append=TRUE )
  
  for ( i in 1:nrow(ctlPars) )
  {
    cat( file=parFile, ctlPars[i,"parameter"]," ",ctlPars[i,"value"],"\n",
         sep="", append=TRUE )
  }

  return( invisible() )
}     # END function saveSimPars

# createList   (converts a data frame into a possiby nested list object)
# Purpose:      Converts a dataframe containing a column "parameter" and a
#               column "value" into a (possibly) nested list.
# Parameters:   obj is a dataframe with columns "parameter" and "value".
# Returns:      result, a list with (possibly) nested levels.
# Source:       A.R. Kronlund
createList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as parameter and values in "value".

  result <- list()

  # Shut off whining, then coerce to numeric to let NA indicate non-numerics.
  options( warn=-1 )
  numericVal <- as.numeric( obj[,"value"] )
  options( warn=0 )

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
      listText <- paste( "result$",obj[i,"parameter"],"=",
                    obj[i,"value"],sep="" )
    # Value is character, build the parse string.
    else
      listText <- paste( "result$",obj[i,"parameter"],"=",
                  obj[i,"value"], sep="" )

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )

    # Evaluate the parse string.
    eval( parse( text=listText ) )
  }
  result
}     # function .createList


# readParFile   (reads an ASCII file with 1 comment line, header, data frame)
# Purpose:      Reads an ASCII file: 1 comment, 1 header, space-delimited
#               data frame usually containing columns "parameter" and "value".
# Parameters:   parFile is a character string indicating the input file.
# Returns:      result, a data frame.
# Source:       A.R. Kronlund
readParFile <- function( parFile="inputFile.par" )
{
  # Read the file and store as a dataframe.

  result <- read.table( file=parFile, as.is=TRUE, header=TRUE, comment.char="#",
                        quote="",sep=" " )
  result
}     # function readPar File



# ctlCreate()
# Function to create a mseR simCtl file from given rep and dat files and
# a base control file.
# inputs:     rep = list created by lisread of a report file
#             age = list created by lisread of ages dat file
#             baseCtl = base control file
#             outCtl = name of txt file for ouput control file.
# outputs:    outCtl = a list object containing the control parameters, also
#                       written to a textfile in the WD.
# usage:      create a new control file for mseR from ADMB output
# source:     SDN Johnson.
ctlCreate <- function ( rep = repList, age = ageList, 
                        baseCtl = "simCtlFile.txt",
                        outCtl = "simCtlFile.txt",
                        simHist = "'simHistNCod.csv'",
                        yrs = tExt, fyr = NULL, proj = 20 )
{
  # Read in base control parameters and coerce to a nested list
  parList   <<- readParFile ( parFile = baseCtl )

  # Compute the time extent from tExt
  if ( !is.null(fyr) ) yrs[1] <- fyr
  tHist <- yrs [ 2 ] - yrs [ 1 ] + 1

  # Now start updating elements of ctlList
  # A function to do this quickly
  replaceParam <- function ( param, new )
  {
    dfRow <- which ( parList [ , "parameter" ] == param )
    parList [ dfRow, "value" ] <<- new
  }

  # First opMod changes
  replaceParam ( param = "opMod$B0", new = rep $ SSB0 / 1000 )
  replaceParam ( param = "opMod$tMP", new = tHist + 1 )   
  replaceParam ( param = "opMod$nT", new = tHist + proj )
  replaceParam ( param = "opMod$historyType", new = simHist )
  # Growth parameters
  replaceParam ( param = "opMod$vonK", new = 0.25 )
  replaceParam ( param = "opMod$Linf", new = 74 )
  replaceParam ( param = "opMod$c1", new = c1 )  # Adjust for mseR units
  replaceParam ( param = "opMod$c2", new = c2 )

  # Age structure
  nAges <- rep $ plusGroupAge
  replaceParam ( param = "opMod$nAges", new = nAges )
  
  # Maturity
  # First, compute maturity schedule
  mat <- rep $ matAge [ , ncol ( rep $ matAge ) ]
  matSched <- approx ( x = mat, y = 1:nAges, xout = c(0.5,0.95) )
  # Put into ctlList
  replaceParam ( param = "opMod$aMat50", new = matSched $ y [ 1 ] )
  replaceParam ( param = "opMod$aMat95", new = matSched $ y [ 2 ] )

  # Selectivity
  replaceParam ( param = "opMod$aSel50", new = rep $ S50_15 )
  replaceParam ( param = "opMod$aSel95", new = rep $ S95_15 )
  replaceParam ( param = "opMod$aSelS50", new = rep $ S50_21 )
  replaceParam ( param = "opMod$aSelS95", new = rep $ S95_21 )

  # Steepness
  replaceParam ( par = "opMod$rSteepness", new = rep $ h )

  # Mortality
  replaceParam ( par = "opMod$M", new = 0.2 )

  # Management procedure parameters which are dependent on time:
  # Assessment
  replaceParam ( param = "mp$assess$t1Method", new =  tHist + 1 )
  replaceParam ( param = "mp$assess$t2Method", new =  tHist + proj )

  # data
  replaceParam ( param = "mp$data$t2Ages", new =  tHist + 1 )
  replaceParam ( param = "mp$data$t2AgesS", new =  tHist + 1 )
  replaceParam ( param = "mp$data$t2Survey", new =  tHist + 1 )

  # hcr
  replaceParam ( param = "mp$hcr$FtEndYr", new =  tHist )
  replaceParam ( param = "mp$hcr$lowBaseEndYr", new =  tHist )
  replaceParam ( param = "mp$hcr$upperBaseEndYr", new =  tHist )

  saveSimPars ( ctlPars = parList, overWrite = TRUE, parFile = outCtl )

  parList
}


# assessCAcreate()
# Function to produce assessCA*.dat files for mseR's assessCA ADMB
# assessment estimator.
# inputs:   repList = list resulting from lisread of rep file
#           ageList = list resulting from lisread of ages dat file
#           idxList = list resulting from lisread of index dat file
#           katList = list resulting from lisread of catch dat file
#           baseCA = list resulting from lisread of base assessCA.dat
#           outFile = character vector name of ouput file
# outputs:  assessCA = list to write out to assessCA*.dat
# usage:    writes a new assessCA*.dat file for use in assessCA ADMB
#           estimator
# source:   SDN Johnson
assessCAcreate <- function (  age, idx, kat, baseCA,
                              outFile = "assessCAout.dat", parList,
                              yrs = c(1959,2014), fyr = 1959 )
{
  # Get years index
  years <- yrs[1]:yrs[2]

  # Compute indices of matrices for time extent
  tIdx <- which ( years == fyr ):length( years )

  # coerce the given parList from data.frame to a list object
  parList <- createList ( parList )

  # Now start the output file
  cat( file=outFile, "## assessCA model data.\n" )
  cat( file=outFile, "## Written:", date(),                     "\n", append=T )
  
  # Now start updating the baseCA list
  # Not all lines need to be updated. Just time series length,
  # and growth parameters to start with
  # First line: time series length and age structure
  baseCA [[ 1 ]] <- c ( parList$opMod$tMP-1,       # Time periods
                        parList$opMod$nAges,    # number of ages
                        1,                      # projection years
                        5,                      # trend years
                        123,                    # rseed
                        2 )                     # hcrType

  # Now growth/maturity parameters
  baseCA [[ 2 ]] <- c ( parList$opMod$aMat50,
                        parList$opMod$aMat95,
                        parList$opMod$Linf,
                        parList$opMod$L1,
                        parList$opMod$vonK,
                        parList$opMod$c1,
                        parList$opMod$c2,
                        parList$opMod$sigmaL )

  # Next catch time series
  baseCA [[ 14 ]] <- katList [[ 1 ]] [ ,3 ] [ tIdx ]

  # Now survey index series
  baseCA [[ 16 ]] <- idxList [[ 8 ]] [ tIdx ]

  # Fishery Age Proportions
  baseCA [[ 17 ]] <- t ( ageList$series1 ) [ tIdx, ]
  # Survey Age Proportions
  baseCA [[ 18 ]] <- t ( ageList$series2 ) [ tIdx, ]

  # Write each element of x to ADMB control file
  cat ( paste ( "Writing to ", outFile, "\n", sep = "") )
  lapply( X=seq_along(baseCA), FUN=writeADMB, x=baseCA, activeFile=outFile )

  baseCA
}


##### Call to create simHistory and simCtl Files #####

# read in repfile and ages file
ageList <- lisread ( agefile )
repList <- lisread ( repfile )
idxList <- lisread ( idxfile )
katList <- lisread ( katfile )
baseCA <- lisread ( "assessCA.dat" )

# Write simHistory and simCtl files to working directory and save a copy
# in a local variable
# 1959
simHist59 <- histCreate ( ages = ageList, rep = repList,
                          yrs = tExt, fyr = 1959,
                          outFile = "simHistNCod.csv" )

simCtl59 <- ctlCreate ( rep = repList, age = ageList, baseCtl = ctlFile,
                        outCtl = "simCtlNCod59.txt", yrs = tExt )

assessCA59 <- assessCAcreate (  age = ageList, baseCA = baseCA,
                                kat = katList, outFile = "assessCA59.dat",
                                parList = simCtl59, fyr = 1959 )

# # 1983
simHist83 <- histCreate ( ages = ageList, rep = repList,
                          yrs = tExt, fyr = 1983,
                          outFile = "simHistNCod83.csv" )

simCtl83 <- ctlCreate ( rep = repList, age = ageList, baseCtl = ctlFile,
                        outCtl = "simCtlNCod83.txt", yrs = tExt,
                        fyr = 1983, simHist = "'simHistNCod83.csv'" )

assessCA83 <- assessCAcreate (  age = ageList, baseCA = baseCA,
                                kat = katList, outFile = "assessCA83.dat",
                                parList = simCtl83, fyr = 1983 )
 


# # 3. Now we want to extract the selectivity at age and maturity
# # schedule info. There's probably more that needs to go into a
# # simCtlFile from the rep file, too
# #   ToDo:   1. Functionalise
# #           2. Read through simCtlFile for ideas on what can come out.
# # Okay, let's pull out the maturity schedule and spline it.
# # matAge <- rep $ matAge

# # matAge <- apply ( X = matAge, FUN = mean, MARGIN = 1 )

# # # Now, use approx to pull out age at 50 and 95 percent mature
# # matSched <- approx ( x = matAge, y = 1:length(matAge), xout = c(0.5, 0.95) )
