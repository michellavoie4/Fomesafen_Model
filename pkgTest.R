###########################################################################################
# A function checking if packages are already installed and if not, install them.

######################################################################################

# INPUT : x : a character vector of package names with 1 or more elements.

# OUTPUT: print if packages are installed. If not, this install the packages.
# If the package name does not exist, a warning is printed : package is not available

pkgTest <- function(x) {
  if(is.null(x)) 
    stop("package names must be inputted")
  if ( is.numeric(x) )
    stop("arguments should be characters")   
  i = 0
  for (i in 1:length(x)) {
    cat('iteration number:', i,'\n')
    cat('package name: ', x[i], '\n')
    if (x[i] %in% installed.packages()) {
      print("Package is already installed") 
      next
    }
    else if (!x[i] %in% installed.packages()) {
      install.packages(x[i],dep=TRUE)
    }
  }
}

