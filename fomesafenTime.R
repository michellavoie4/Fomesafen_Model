###########################################################################################
# A script calculating the time for code execution
############################################################################################

nb <- 20                  # number of time the code is timed 
user_time <- numeric(nb)  # Initialize the vector of total user time for code profiling
k <- 0
for (k in 1:nb) {
  user_time[k] <- system.time({
    source(paste0(getwd(),"/mainFomesafen.R")) })[1] # Execute the code and store user time
}

message("--------------------------------------------------------")
message("########### Code Benchmarking ############")
cat("Number of time the code was executed : ", nb, '\n')
cat("Mean elapsed time is : ", mean(user_time), "s", '\n')
cat("Standard deviation of elapsed time is : ", sd(user_time), "s", '\n')
