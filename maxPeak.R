###########################################################################################################################
# A function calculating : 1) the maximum cellular concentration of fomesafen (mol/m^3) (maxCell)
# 2) Time (h) needed to reach maximum cellular fomesafen concentration (tPeak)
# and summarizing the key results on fomesafen uptake kinetics
###########################################################################################################################

# INPUT : A matrix giving cellular fomesafen concentrations for the 4 different fomesafen concentrations (in columns) over time (first row)

# OUTPUT: This function returns a dataframe with total fomesafen concentration, protonated fomesafen concentration,
#           maximum cellular fomesafen concentration at low or high cell radius, time needed to reach maximum cellular fomesafen concentration

maxPeak <- function(x) {
  if( !is.matrix(out) || !dim(out)[2] == 10 ) 
    stop("The function argument should be a matrix with 10 columns")
  
  i <- 0
  maxCell <- c(0)
  tPeak <- c(0)
  for (i in 1:4) {
    maxCell[i] = max(x[,i+1])
    tPeak[i] = dt*(which.max(x[,i+1]) - dt)
  }
  maxCell
  tPeak
  
  i <- 0
  maxCell_high <- c(0)
  tPeak_high <- c(0)
  for (i in 1:4) {
    maxCell_high[i] = max(x[,i+6])
    tPeak_high[i] = dt*(which.max(x[,i+6]) - dt)
  }
  maxCell_high
  tPeak_high
  rmax <- maxCell/maxCell_high
  rt <- tPeak_high/tPeak
  
  return(data.frame("Total_fomesafen_Conc" = Conc_ug, "Protonated_Fomesafen_Conc" = Conc_AH, "max_Conc_low_R" = maxCell, "max_Conc_high_R" = maxCell_high, "Conc_Ratio" = rmax, "Time_Peak_lowR" = tPeak,  "Time_Peak_highR" = tPeak_high, "Time_Ratio" = rt))
}
###########################################################################################################################
