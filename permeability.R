####################################################################################################################
# Function calculating cell permeability and protonated fraction of fomesafen
####################################################################################################################

# INPUT : pH : pH of the medium
#         pKa : pKa of the toxicant (fomesafen)
#         R : algal cell radius (in meter)

# OUTPUT: return a list Papp : apparent cellular permeability (m/h)
#                       Conc_AH : protonated fomesafen concentration at pH 6.8 (mol / m^3)
#                       Pm and Pbl : membrane and boundary layer permeability (m/h)
#                       Pm_Pbl_ratio : Ratio of Pm to Pbl

permeability <- function(pH, pKa, R) {
  if(is.null(pH) || is.null(pKa) || is.null(R)) 
    stop("all arguments must be inputted")
  if ( !is.numeric(pH) || (!is.numeric(pKa)) || (!is.numeric(R)) )
    stop("all arguments should be numeric") 
  if ( pH < 0 || pH > 14 || pKa < 0 || pKa > 14 )
    stop("pH and pKa should be between 0 and 14") 
  
  # Computing the concentration of the protonated fomesafen species using the Henderson-Hasselbalch equation
  A_AH_ratio <- 10^(pH - pKa)            # Ratio of unprotonated to protonated species 
  F_AH <- 1-(A_AH_ratio/(1+A_AH_ratio))  # Fraction of protonated molecule (no units)
  Prop_AH <- F_AH * 100                  # proportion of protonated molecule (%)
  Conc_AH <- Conc*F_AH                   # Initial concentration of protonated fomesafen in the culture medium (mol/m^3)
  Conc_AH
  
  # Calculation of fomesafen permeability
  Kow_exp <- 2.90                         # Octanol water partition coefficient from Experimental database match (Tomlin, 1994) 
  # Kow_model <- 3.41                       # Octanol water partition coefficient (KOWWIN v1.67 estimates)
  Pm_cm <- 10^1.11*log10(Kow_exp) - 0.6            # (cm/s) Empiricial equation of lipid bilayer permeability of non-electrolyte : log Pm = s log (kow) + b, where s = 1.11 and b = -0.67 for the solutes with MW > 50 Da (Walter and Gutknecht, 1986)
  Pm <- Pm_cm /100                        # (m/s)
  
  # Calculating the diffusion coefficient of fomesafen
  # kb <- 1.38064852E-23                              # Boltzmann's constant (J/K) (or Kg m^2 s-2 K-1)
  # T <- 25 + 273                                     # Temperature (K) 
  # n <- 0.89                                         # Absolute or dynamic viscosity of water at 25 °•C (N s / m^2) (or: Kg m s-2 s m-2 = Kg s-1 m-1)
  # R_fomesafen <- 1E-12                              # Hydrodynamic radii (m)
  # D <- (kb * T) / (6 * pi * n * R_fomesafen)        # Diffusion coefficient of fomesafen (m^2/s)
  # D
  
  # Calculating the diffusion coefficient of fomesafen with the empirical equation of Hayduk and Laudie (1974)
  cP <- 0.8937                    # viscosity of water at 25 C (centipoise) (Haynes et al., 2016)
  Vm <- 278.7                     # molar volume (cm^3/mol) (LookChem website)
  Dcm <- (13.3E-05) / ((cP ^ -1.14) * (Vm ^ 0.589))  # Diffusion coefficient in cm^2/s (Hayduk and Laudie, 1974)
  D <- Dcm/1E+04                  # Diffusion coefficient (m^2/s)
  
  # Caclulating the the overall permeability of fomesafen in the boundary layer and the cell membrane 
  R                                    # Cell radius
  L <- R                               # Boundary layer thickness (m) equals to cell radius (Wolf-Gladrow and Riebesell, 1997; Lavoie et al, 2018)
  Pbl <- D/L                           # Permeability of the boundary layer (m s-1)
  Papp <- 1/(1/Pm + 1/Pbl)             # Apparent permeability of the membrane plus that of the boundary layer (m s-1) (Missner and Pohl, 2009)
  Pm_Pbl_ratio <- Pm/Pbl               # Ratio Pm/Pbl, ratio of membrane permeability to boundary layer permeability
  
  # Converting the seconds in hours for all relevant parameters
  D1 <- D * 60 * 60                 # D (m^2/h)
  Pbl_1 <- D1/L                      # Permeability of BL (m/h)
  Pm1 <- Pm * 60 * 60                # Permeability of membrane (m/h)
  Papp1 <- (1/(1/Pm1 + 1/Pbl_1))       # Papp (in m/h)
  
  return(list(Papp1=Papp1, Conc_AH=Conc_AH, Pm=Pm, Pbl=Pbl, Pm_Pbl_ratio=Pm_Pbl_ratio))
}
##################################################################################################
