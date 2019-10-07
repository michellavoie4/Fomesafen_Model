# Procedure to run the model predicting cellular fomesafen uptake over time in a phytoplanktonic cell, which is Raphidocelis subcapitata in this case. We consider four different initial dissolved total fomesafen concentrations.

1) Download and install the latest versions of R at https://cran.rstudio.com/

2) Download and install the latest version of RStudio at https://www.rstudio.com/products/rstudio/download/

3) Copy the GitHub repository using the "git clone" command in the terminal (Mac or Linux) or in git bash (Windows) or by downloading it from GitHub.

4) In RStudio, set the new directory called "fomesafen" as your R working directory.

5) Run (i.e., source) the "mainFomesafen.R" file in RStudio, which produces the plots of cellular fomesafen concentrations as a function of exposure time.

6) The "mainFomesafen.R" file uses 4 functions, which are in the same directory than "mainFomesafen.R" and are called "pkgTest.R", "permeability.R", "plotRes.R", and "maxPeak.R".

7) Note that the script "mainFomesafen.R" first check if the required packages (i.e., tidy verse, deSolve, and testit) are already installed on your computer. If they are not, they will be automatically installed  using your default CRAN mirror. If you do not have a default CRAN mirror, R will ask you to select one.

8) Note also that the R script (fomesafenTime.R), which computes the mean time needed to run the model, can also be run (i.e., source) independently for code benchmarking or profiling.
