###########################################################################################################################
# Function to plot the results in 4 panels (1: 5 ug/L, 2: 10 ug/L; 3: 40 ug/L; 4: 320 ug/L)
###########################################################################################################################

# INPUT : A matrix giving cellular fomesafen concentrations for the 4 different fomesafen concentrations (in columns) over time (first row)

# OUTPUT : Two plots of the results with 4 panels each. One plot assuming no fomesafen loss and the other assuming passive eflux of fomesafen

plotRes <- function(out) {
  if( !is.matrix(out) || !dim(out)[2] == 10 ) 
    stop("The function argument should be a matrix with 10 columns")
  
  out <- out[which(out[,1] < 1.1),] # Only select the first 1.1 minute
  times <- times[which(times < 1.1)] # Only select the first 1.1 minute
  
  old.par <- par()
  mf <- par(mfrow = c(2, 2))
  matplot(out[,1], out[,2], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,2])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  lines(out[,1], out[,7], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,2]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
  text(0.1, max(out[,2])*1.1, "A", cex=1.5)
  matplot(out[,1], out[,3], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,3])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
  lines(out[,1], out[,8], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,3]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  text(0.1, max(out[,3]*1.1), "B", cex=1.5)
  matplot(out[,1], out[,4], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,4])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
  lines(out[,1], out[,9], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,4]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  text(0.1, max(out[,4])*1.1, "C", cex=1.5)
  matplot(out[,1], out[,5], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,5])*1.2), type="l", xlab="Time (min)", ylab="Cell concentration (mol/m^3)")
  lines(out[,1], out[,10], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,5]*1.2)), xlab=c("Time (min)"), ylab="cell concentration (mol/m^3)")
  axis(1,at=c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  text(0.1, max(out[,5])*1.1, "D", cex=1.5)
  par(mf)
  old.par
  #dev.off()
  
  # setEPS()
  # for (k in 1:length(outlist)) {
  #   if (k == 1) {
  #     postscript("out.eps")
  #   } else 
  #     postscript("out_loss.eps")
  # }
  # old.par <- par()
  # mf <- par(mfrow = c(2, 2))
  # i <- 0
  # for (i in 1:4) {
  #   matplot(out[,1], out[,1+i], xaxt="n", xlim=c(0,max(times)), ylim=c(0,max(out[,1+i])*1.2), type="l", xlab="Time (h)", ylab="Cell concentration (mol/m^3)")
  #   axis(1,at=c(0,10,20,30,40,50,60,70))
  #   text(5, max(out[,i+1])*1.1, LETTERS[i], cex=1.5)
  #   lines(out[,1], out[,i+6], xlim=c(0,times[length(times)]), lty = 2, ylim=c(0, max(out[,i+1]*1.2)), xlab=c("Time (h)","Time (h)","Time (h)", "Time (h)"), ylab="cell concentration (mol/m^3)")
  # }
  # 
  # par(mf)
  # old.par
  # dev.off
  
}
###########################################################################################################################
