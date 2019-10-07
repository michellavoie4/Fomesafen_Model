###########################################################################################
# A script calculating the time for code execution
############################################################################################

nb <- 2                  # number of time the code is timed 
start_time <- numeric(nb) # initialize the vector of start time for code profiling
user_time <- numeric(nb)  # Initialize the vector of total user time for code profiling
k <- 0
for (k in 1:nb) {
  start_time[k] <- proc.time()[1]  # Store the initial time for each code profiling test
  source(paste0(getwd(),"/fomesafen2.R"))  # Execute the code
  user_time[k] <- (proc.time()[1] - start_time[k][1])  # Stop the clock and print elapsed time
}

message("--------------------------------------------------------")#, fill=T)
message("########### Code Benchmarking ############")#, fill=T)
cat("Number of time the code was executed : ", nb, fill=T)
cat("Mean elapsed time or wall clock time is : ", mean(user_time), fill=T)
cat("SD of elapsed time or wall clock time is : ", sd(user_time), fill=T)

