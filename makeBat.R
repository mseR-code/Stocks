# Little script to create a batch file for LOTS of different MPs in mseR

# Set up experiments to run batches over by creating lists
# of parameters and different values at which to test them.
# Function below also requires abbreviations of the 

# First, let's look at performance of DCAC over different
# multiples of actual depletion and natural mortality, and
# diffeerent levels of uncertainty in those estimates.
parListDM <- list ( 'mp$assess$dcacDeltaMult' = seq ( 0.5, 2, by = 0.5 ),
                    'mp$assess$dcacSDdelta' = seq ( 0.1, 0.5, by = 0.1),
                    'mp$assess$dcacMortMult' = seq ( 0.5, 2, by = 0.5),
                    'mp$assess$dcacCVLogM' = seq ( 0.2, 1, by = 0.2),
                    'mp$assess$dcacQuant' = c(TRUE, FALSE) )
labelDM <- c ( "Dm", "Dv", "Mm", "Mv", "Q" )

parListDMs <- list (  'mp$assess$dcacDeltaMult' = seq ( 0.5, 2, by = 0.5 ),
                      'mp$assess$dcacSDdelta' = seq ( 0.1, 0.5, by = 0.1),
                      'mp$assess$dcacMortMult' = seq ( 0.5, 2, by = 0.5),
                      'mp$assess$dcacCVLogM' = seq ( 0.2, 1, by = 0.2),
                      'mp$assess$dcacQuant' = c(TRUE, FALSE),
                      'mp$assess$dcacState' = TRUE )
labelDMs <- c ( "Dm", "Dv", "Mm", "Mv", "Q", "S" )

# Now let's look at HCR settings - this should be tested
# against different settings above too, but those should be 
# identified for efficient experiments
parListHCRq <- list(  'mp$assess$dcacCons' = seq ( 0.2, 1, by = 0.2 ),
                      'mp$assess$dcacPosGrad' = seq ( 0, 3, by = 0.5 ),
                      'mp$assess$dcacNegGrad' = seq ( 0, 3, by = 0.5 ),
                      'mp$assess$dcacState' = c(TRUE, FALSE),
                      'mp$assess$dcacQuant' = TRUE )
parListHCRc <- list(  'mp$assess$dcacCons' = seq ( 0.2, 1, by = 0.2 ),
                      'mp$assess$dcacPosGrad' = seq ( 0, 10, by = 2 ),
                      'mp$assess$dcacNegGrad' = seq ( 0, 10, by = 2 ),
                      'mp$assess$dcacState' = c(TRUE, FALSE),
                      'mp$assess$dcacQuant' = FALSE  )
labelHCR <- c ( "Pr", "PG", "NG", "St", "Q" )


# And trend detection settings.
parListTrend <- list (  'mp$assess$dcacPosTrendProb' = seq ( 0.1, 0.5, by = 0.2),
                        'mp$assess$dcacNegTrendProb' = seq ( 0.1, 0.5, by = 0.2),
                        'mp$assess$dcacTimeFrame' = seq ( 3, 7, by = 2),
                        'mp$assess$dcacState' = c(TRUE, FALSE),
                        'mp$assess$dcacQuant' = c(TRUE, FALSE) )
labelTrend <- c( "PT", "NT", "TF", "St", "Q" )

# Splitting the experiment up this way cuts down on the combinatorial
# nature of the experiments: doing all combinations at once would have
# led to 84m trials.


# Function to create batch files from the lists above
batCreate <- function ( parList, outFile, label )
{
  exp <- expand.grid ( parList )
  for ( j in 1:nrow ( exp ))
  {
    # Create a label
    lab <- character()
    for ( k in 1:ncol ( exp ) ) 
      lab <- paste ( lab, label[k], as.numeric(exp[j,k]), sep = "" )
    # Print MP name
    if ( j == 1 ) 
      cat ( "# Management Procedure ", j, " : ", lab,
            "\n", sep = "", file = outFile, append = FALSE)
      else cat (  "# Management Procedure ", j, " : ", lab, 
                  "\n", sep = "", file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
      cat ( "mp$mp", j, "$gui$mpLabel '", lab, "'\n", sep = "", append = TRUE, 
            file = outFile )
      for ( k in 1:ncol ( exp ) )
        cat ( "mp$mp", j, "$", names(exp)[k], " ", exp[j,k], "\n", sep ="",
              file = outFile, append = TRUE )
      cat ( "#\n", file = outFile, append = TRUE )
  }
  cat ( "# File Ends <not run>.", file = outFile, append = TRUE )
}



