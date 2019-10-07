###############################################################################################
# Solving the differential equation expressing the rate of change in cellular
# fomesafen over time for 4 different dissolved fomesafen concentrations
#
# Author: Michel Lavoie
#
# Last update: October 6 2019
#
# This code has been written and tested with R version 3.6.0 (version$version.string)
# I used packages "deSolve' v 1.21, "testit" v 0.9, and "tidyverse" v 1.2.1 (packageVersion())
###############################################################################################

# First set your working directory to the folder 'Fomesafen'

###############################################################################################
# Check for installed packages and install them if they are not installed.
pkg <- c("deSolve", "testit", "tidyverse")
source(paste0(getwd(),"/pkgTest.R"))
pkgTest(pkg)

# loading packages
lapply(pkg, library, character.only = T)
###############################################################################################


###############################################################################################
# parameters and initial calculations
MW <- 438.76                              # Molecular weight of fomesafen (g/mol) (Pubchem website)
Conc_ug <- c(5, 10, 40, 320)              # Concentration in ug/L
Conc_molL <- (Conc_ug/MW) / 1E+06         # Concentration in mol/L
Conc <- Conc_molL*1000                    # Concentration in mol/m^3
Sol_mg <- 50                           # (mg/L) Water solubility measured at 20 C (Shiu et al., 1990)
Sol <- Sol_mg/1000/MW*1000               # (mol/m^3) Water solubility
assert(max(Conc) < Sol)                  # Make sure all fomesafen concentrations are lower than the solubility limit
cat("Ratios of fomesafen concentrations to solubility limit :", Conc/Sol, '\n')
pKa <- 3.0                                # pKa (Weber, 1993)
pH <- 6.8                                  # Initial measured pH of the culture medium (This study)
assert(Conc & Sol & pKa & pH > 0)         # Make sure all parameters are positive (no sign errors)

# Validating that the proportion of culture medium volume occupied by the algae is small
cell_dens <- 1000000             # Macimum final cell density (cell/mL)
biov_cell <- 34*1E-15            # Maximum cell biovolume (m^3 / cell)
biovL <- cell_dens * biov_cell * 1000  # total biovolume (m^3 / L of medium)
biovm <- biovL * 1000                    # total biovolume (m^3 / m^3 of medium)
biovL <- biovL * 1000                # total biovolume (L / L of medium) 1 cubic meter = 1000 liters
cat("Naximum proportion of the culture medium volume occupied by the algae : ", biovL * 100, "%", '\n')

# Validating that the proportion of culture medium volume occupied by the algae is small
cell_dens <- 116280             # Mean measured initial cell density (cell/mL) (n=45)
biov_cell <- 30.33*1E-15            # Mean measured initial cell biovolume (m^3 / cell) (n=45)
biovL <- cell_dens * biov_cell * 1000  # total biovolume (m^3 / L of medium)
frac <- biovL * 1000                    # total biovolume (m^3 / m^3 of medium)  # Mean measured initial volume fraction of algae in the mediu
cat("Initial fraction of the culture medium volume occupied by the algae : ", frac, "(no units)", '\n')
###############################################################################################


############################################################################################################################
# Variables
u1 <- c(0.72, 0.74, 0.71, 0.51)                 # Mean measured growth rate of R. subcapitata (d-1) (This study)
u1 <- u1 / 24                                    # Mean measured growth rate of R. subcapitata (h-1) (This study)
Ri <- c(1.92,1.93,1.93,1.92)                         # Mean measured initial cell radius of R. subcapitata (um) (This study)
Ri <- Ri * 1E-06                                     # Mean measusred initial cell radius (m)
Rf <- c(2.02,2.05,2.21,2.33)                                   # Mean measured FINAL cell radius of R. subcapitata (um) (This study)
Rf <- Rf * 1E-06                                      # Mean measusred FINAL cell radius (m)
Vi <- 4/3 * pi * (Ri^3)                 
Ai <- 4 * pi * (Ri^2)
Aini <- Ai / Vi                                   # Mean measured INITIAL cellular surface : volume ratio (m^2/m^3)
Vf <- 4/3 * pi * (Rf^3)                 
Af <- 4 * pi * (Rf^2)
Afi <- Af / Vf                                    # Mean measured FINAL cellular surface : volume ratio (m^2/m^3)
A_R <- rbind(Aini, Afi)                           # Preparing the matrix of cellular surface : volume ratio
############################################################################################################################


############################################################################################################################
# Calculations of cell permeability and protonated fomesafen fraction
# This is needed for the derivative equation below
############################################################################################################################
source(paste0(getwd(),"/permeability.R"))
Pbl_R <- rbind(permeability(pH=pH, pKa=pKa, R=Ri)$Pbl, permeability(pH=pH, pKa=pKa, R=Rf)$Pbl)   # Cellular permeability considering a low cell radius (1st row) or a high cell radius (2sd row)
Conc_AH <- permeability(pH=pH, pKa=pKa, R=Ri)$Conc_AH                     # Does not change with cell radius
assert(u1 & Ri & Rf & frac & as.vector(Pbl_R) & Conc_AH > 0) # Papp_R & Conc_AH > 0)                 # Make sure all variables are positive (no sign errors)

cat("Calculated membrane permeability (Pm) in m/h at low or high cell radius : ", permeability(pH=pH, pKa=pKa, R=Ri)$Pm, fill=T) # Equal at low or high cell radius
cat("Ratio of membrane permeability to boundary layer permeability at low cell radius:", permeability(pH=pH, pKa=pKa, R=Ri)$Pm_Pbl_ratio, fill=T)
cat("Ratio of membrane permeability to boundary layer permeability at high cell radius:", permeability(pH=pH, pKa=pKa, R=Rf)$Pm_Pbl_ratio, fill=T)
############################################################################################################################



##############################################################################################################
# Solving the differential equations expressing the rate of change of cellular fmesafen over time
#############################################################################################################

t <- 1           # exposure time (in hours)
dt <- 0.00003# 0.0003       # time step for numerical integration (in h)
cat("time step (in sec) : ", dt*60*60, '\n')
cat("Number of time steps : ", t/dt, '\n')
out <- matrix(0, nrow = t/dt+1, ncol = 0)    # initialize the matrix of cellular fomesafen concentrations assuming no efflux
out_loss <- out                              # initialize the matrix of cellular fomesafen concentrations assuming efflux
outnew <- matrix(0, nrow = t/dt+1, ncol = 0) # initialize the temprory matrix of cellular fomesafen concentrations assuming no efflux
outnew_loss <- outnew                        # initialize the temporay matrix of cellular fomesafen concentrations assuming efflux

i <- 0
for (i in 1:2) {

# Differential equation expressing the rate of change of cellular fomesafen concentration (dy, units = mol / m^3)
# function using low cell radius (i = 1) and high cell radius (i = 2)
cell_fomesafen <- function(t, y, parms) {
  dy <- (A_R[i,] * Pbl_R[i,] * ( Conc_AH - y) ) - (u1 * y) 
  #dy <- (A_R[i,] * Papp_R[i,] * (Conc_AH - (y * ( frac * exp(1)^(u1 * t)) ) - y )) - (u1 * y)  # Considering depletion of Conc_AH in the medium
  list(dy)
}

# Differential equation taking into account active fomesafen excretion or metabolism
cell_fomesafen_loss <- function(t, y, parms) {
  ke <- A_R[i,] * Pbl_R[i,]
  dy <- (A_R[i,] * Pbl_R[i,] * (Conc_AH - y) )  - ( ke * y ) # - (u1 * y)
  # dy <- (A_R[i,] * Pbl_R[i,] * (Conc_AH - y) )  - ( A_R[i,] * - Pbl_R[i,]  * (y - Conc_AH) ) # Gives zero.. assuming passive efflux
  list(dy)
}

# Solving dynamically the cellular fomesafen concentration as a function of 4 external fomesafen concentrations
yini <- c(0,0,0,0)  # initial dissolved fomesafen concentration
times <- seq(from = 0, to = t, by = dt)  
outnew <- ode (times = times, y = yini, func = cell_fomesafen)
outnew_loss <- ode (times = times, y = yini, func = cell_fomesafen_loss)
out <- cbind(out, outnew)  # Matrix : in rows = time points; in columns : 1: elapsed time, 2: 5 ug/L, 3: 10 ug/L, 4: 40 ug/L, 5: 320 ug/L.
out_loss <- cbind(out_loss, outnew_loss) 
}

out[,1] <- out[,1] * 60               # convert time in the matrix 'out' in minutes
out_loss[,1] <- out_loss[,1] * 60
times <- times * 60                    # Convert the variable "times" in minutes
##########################################################################################################################


###########################################################################################################################
# Printing a data frame summarizing the results : maximum fomesafen concentration and time needed to reach maximum for both set of plots assuming low or high cell radius.
# Case 1 : Without fomesafen efflux
# Case 2: With fomesafen efflux
outlist <- list(out, out_loss)
names(outlist) <- c("Case 1", "Case 2")
source(paste0(getwd(), "/maxPeak.R"))
map(outlist, maxPeak) # Calling the function "max_Peak"
###########################################################################################################################


###########################################################################################################################
# Plot the results for the case when fomesafen efflux is neglected : 5, 10, 40, 320 ug/L fomesafen
setEPS()
postscript("out.eps")
source(paste0(getwd(), "/plotRes.R"))
plotRes(out) # call the function "plotRes"
dev.off()

# Plot the results when fomesafen efflux is considered : 5, 10, 40, 320 ug/L fomesafen
setEPS()
postscript("out_loss.eps")
plotRes(out_loss)
dev.off()

###########################################################################################################################


######################################################################################################
# Calculating the effect of ke on fomesafen uptake kinetics
# Differential equation taking into account active fomesafen excretion or metabolism
# Using a fix low initial cell radius

t <- 1           # exposure time (in hours)
dt <- 0.00003# 0.0003       # time step for numerical integration (in h)
cat("time step (in sec) : ", dt*60*60, '\n')
cat("Number of time steps : ", t/dt, '\n')
out <- matrix(0, nrow = t/dt+1, ncol = 0)    # initialize the matrix of cellular fomesafen concentrations assuming no efflux
ke <- c(0.1, 30, 300)   # fomesafen loss rate constant (h-1)
out_loss2 <- matrix(0, nrow = t/dt+1, ncol = 0)    # initialize the matrix of cellular fomesafen concentrations assuming no efflux

i <- 0
for (i in 1:length(ke)) {

cell_fomesafen_loss2 <- function(t, y, parms) {
  dy <- (A_R[1,] * Pbl_R[1,] * (Conc_AH - y) )  - ( ke[i] * y ) # - (u1 * y)
  # dy <- (A_R[i,] * Pbl_R[i,] * (Conc_AH - y) )  - ( A_R[i,] * - Pbl_R[i,]  * (y - Conc_AH) ) # Gives zero.. assuming passive efflux
  list(dy)
}

# Solving dynamically the cellular fomesafen concentration as a function of 4 external fomesafen concentrations
yini <- c(0,0,0,0)  # initial dissolved fomesafen concentration
times <- seq(from = 0, to = t, by = dt)  
outnew_loss2 <- ode (times = times, y = yini, func = cell_fomesafen_loss2)
out_loss2 <- cbind(out_loss2, outnew_loss2) 
}

# Plot the results
out_loss2[,1] <- out_loss2[,1] * 60               # convert time in the matrix 'out' in minutes
times <- times * 60                    # Convert the variable "times" in minutes

out_loss2 <- out_loss2[which(out_loss2[,1] < 1.1),] # Only select the first 1.1 minute
times <- times[which(times < 1.1)] # Only select the first 1.1 minute
out <- out_loss2

setEPS()
postscript("out_ke.eps")

old.par <- par()
mf <- par(mfrow = c(2, 2))
  
matplot(out[,1], out[,2], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,2])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
lines(out[,1], out[,7], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,3]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
lines(out[,1], out[,12], xlim=c(0,times[length(times)]), lty = 3, ylim=c(0, max(out[,3]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
text(0.1, max(out[,2])*1.1, "A", cex=1.5)

matplot(out[,1], out[,3], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,3])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
lines(out[,1], out[,8], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,3]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
lines(out[,1], out[,13], xlim=c(0,times[length(times)]), lty = 3, ylim=c(0, max(out[,3]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
text(0.1, max(out[,3]*1.1), "B", cex=1.5)

matplot(out[,1], out[,4], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,4])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
lines(out[,1], out[,9], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,4]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
lines(out[,1], out[,14], xlim=c(0,times[length(times)]), lty = 3, ylim=c(0, max(out[,4]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
text(0.1, max(out[,4])*1.1, "C", cex=1.5)

matplot(out[,1], out[,5], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,5])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
lines(out[,1], out[,10], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,5]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
lines(out[,1], out[,15], xlim=c(0,times[length(times)]), lty = 3, ylim=c(0, max(out[,5]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
text(0.1, max(out[,5])*1.1, "D", cex=1.5)

par(mf)
old.par

dev.off()


#########################################################################################################
#########################                 References             #######################################
########################################################################################################

# Wolf-Gladrow, D.; Riebesell, U., Diffusion and reactions in the vicinity of plankton: A refined model for inorganic carbon transport. Mar. Chem. 1997, 59 (1-2), 17-34.
# Lavoie, M.; Galí, M.; Sévigny, C.; Kieber, D. J.; Sunda, W. G.; Spiese, C. E.; Maps, F.; Levasseur, M., Modelling dimethylsulfide diffusion in the algal external boundary layer: implications for mutualistic and signalling roles. Environ. Microbiol. 2018, 20 (11), 4157-4169.
# Pubchem website. Page visited on May 22th 2019. https://pubchem.ncbi.nlm.nih.gov/compound/Fomesafen#section=Melting-Point
# Shiu, W.; Ma, K.; Mackay, D.; Seiber, J.; Wauchope, R., Solubilities of pesticide chemicals in water. Part II: Data compilation. Rev. Environ. Contam. Toxicol. 1990, 116: 15-187.
# Weber, J. B., Ionization and sorption of fomesafen and atrazine by soils and soil constituents. Pestic. Sci. 1993, 39 (1), 31-38.
# Tomlin, C. (ed.). 1994. The pesticide manual: A world compendium. 10th ed. (Incorporating the Agrochemicals handbook.) British Crop Protection Council and Royal Society of Chemistry, Thornton Heath, UK
# Walter, A.; Gutknecht, J., Permeability of small nonelectrolytes through lipid bilayer membranes. J. Membrane Biology 1986, 90, 207-217.
# Haynes, W. M.; Lide, D. R.; Bruno, T. J., CRC handbook of chemistry and physics : a ready-reference book of chemical and physical data. 97th Edition ed.; Boca Raton, CRC Press: Florida, 2016.
# LookChem website, page visited on May 22th 2019. https://www.lookchem.com/Fomesafen/
# Missner, A.; Pohl, P., 110 Years of the Meyer-Overton rule: predicting membrane permeability of gases and other small compounds. Chem Phys Chem 2009, 10 (9-10), 1405-1414.
# Hayduk, W.; Laudie, H., Prediction of diffusion coefficients for nonelectrolytes in dilute aqueous solutions. AICHE J. 1974, 20 (3), 611-615.


###########################################################################################################

# NOTES
# At lower cell radius, the peak fomesafen cellular concentration is reached more quickly (after 28.9 h rather than 30.7 h)
# Therefore, increase in cell radius will have only a small protectory effects on fomesafen uptake/toxicity .

# What is the total fomesafen amount taken up after 72 h?
# What is the instantaneous product of cell density and fomesafen quotas at each time point (total fomesafen taken up)?
# What is the mean or range in cellular fomesafen concentrations between 30 and 72 hours? Perhaps a decrease in growth rate at 320 ug/L fomesafen is beneficial and decrease much fomesafen uptake

# Decreasing growth rate, delay reaching maximum cellular fomesafen concentration (by several hours), but also delay depuration.
# Interestingly, in all cases (at pH 6.8) , time was sufficient for reaching high intracellular fomesafen concentrations (comparable to thosein the medium),
# but, an increase in pH will also decrease uptake of fomesafen over time and further delay reaching cellular fomesafen peak

# Postponing cellular fomesafen depuration when the growth rate decrease (as well as cell radius increases) will increase the apparent toxicity after 72 h, but decrease it on the short term (e.g., 24h)
# This trend over time fo fomesafen toxicity could be further amplified with a rise in pH over time (due to photosynthesis).

# So, the fomesafen toxicity measured after 72 h of exposure to 320 ug/L should be proportionatly more important than at 40 ug/L!!
# Do we see this in the data?
# YES !! We see that the mean growth rate decreases 1.8-fold between 40 and 320 ug/l (for a 8-fold increase in fomesafen concentration!)
# But, we see that the mean growth rate decreases only by 20% between 10 and 40 ug/L fomesafen !!!!!

# YES !! We also see that mean maximum and operational PSII activity decreases by 10 and 22% for a 8-fold increases in fomesafen concentration.
# but, we see that mean maximum and operational PSII activity decreases by only 3 and 5% for a 4-fold increases in fomesafen (from 10 to 40 ug/L)
# Therefore, the model demonstrates that negative feedback exist on fomesafen toxicity due to time delay in fomesafen uptake and depuration, which is due to a decrease in growth rate and increase in cell radius.

# Perhaps, I should make 4 plots with respective mean growth rate and mean final cell radius.
# And make 4 other plots with respective growth rate and mean INITIAL cell radius.
# This should illustrate the point above.


# There are no steady-state reached since fomesafen  is only added in the medium at the start of the experiment.
# If we would have use a chemostat, a steady-state would have been reach which would have shown that the steady-state fomesafen concentration
# decrease as we decreases the dissolved fomesafen concentration or decrease cell radius...
