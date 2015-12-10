# --------------------------------------------------------------------
# compareOMRep.R
#
# A script to compare mseR's operating model output to the output
# of the ADMB model assessment's rep file.
# 
# Compares depletion of SSB at tMP-1, and plots the time series of
# Bt from the assessment and from the OM.
#
# --------------------------------------------------------------------

# Source mseRtools
source ("mseRtools.R")

# load rep file
rep <- lisread ( "ncod.rep" )

# recover rep file quantities/series
repSSBt <- rep$SSBt / 1000            # Convert to 1e3t
nT <- length ( repSSBt )              # ts length
repDep <- rep$SSBt [ nT ] / rep$SSB0  # compute depletion

# Load blob (unruly filenames, ugh, probably drag filepath in)
load("/Users/samuelj/Work/code/mseR/mseR-NCod/mseRproject/sim_27112015132039/sim_27112015132039.Rdata")

# Take the number of years used in the operating model history
OMnT <- blob$ctlList$opMod$tMP-1

# Recover Bt series
omBt <- blob$om$Bt[1,2:(OMnT+1)]
omDep <- omBt[OMnT] / blob$ctlList$opMod$B0

# Compute residuals
residBt <- omBt - repSSBt[(nT - OMnT + 1):nT]
residDep <- round ( omDep - repDep, digits = 4)

# Plot states and residuals
par ( mfrow = c(2,1))
plot ( x = c(1,OMnT), y = c ( 0, max ( omBt, repSSBt) ), type = "n" )
  lines ( omBt, lty = 2, lwd = 2, col = "grey70" )
  lines ( repSSBt[(nT - OMnT + 1):nT], lty = 1, lwd = 2 )
  mtext ( text = "SSB")
  panLegend ( legTxt = c("rep", "om"), col = c("black", "grey70"),
              lty = c(1,2), x = 0.8, y = 1, bty = "n" )
  panLab ( x = 0.2, y = 0.8, txt = paste ( "Residual Dep = ", residDep, sep = "" ) )
barplot ( residBt, col = "grey70", main = "Residuals" )

