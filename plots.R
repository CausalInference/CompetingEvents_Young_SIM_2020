###########
###### Plots. Requires all other programs to be run first
###########
setwd("your directory")
pdf("Risk5Event.pdf")
plot(cutTimes,cumIncTreated5, type="s",ylim=c(0,1), ylab="Risk", 
     xlab="Month",lty=1,lwd=1,xlim=c(0,60), cex.lab=1.2,cex.main=1.2,main =c(paste("Risk of prostate cancer death"),paste("without elimination of competing events"))) # Death due to Prostate Cancer: Estimand (5)
lines(cutTimes, cumIncPlacebo5, type="s",col=2,ylim=c(0,1),lwd=1)
lines(cutTimes, cumIncTreatedIpw5_2, type="s",col=1,lwd=1,lty=2)
lines(cutTimes, cumIncPlaceboIpw5_2, type="s",col=2,lwd=1,lty=2)
lines(cutTimes, cumIncTreatedIpw5, type="s",col=1,ylim=c(0,1),lwd=1,lty=3)
lines(cutTimes, cumIncPlaceboIpw5, type="s",col=2,ylim=c(0,1),lwd=1,lty=3)
legend("topleft", c("parametric g-formula (treatment)","parametric g-formula (placebo)","IPW sub (treatment)","IPW sub (placebo)","IPW cs (treatment)","IPW cs (placebo)"),
       col=c(1,2,1,2,1,2), lty=c(1,1,2,2,3,3),cex=1.2,pt.cex=1.2,lwd=2)
dev.off()

pdf("Risk3Event.pdf")
plot(cutTimes,cumIncTreated3, type="s",ylim=c(0,1), ylab="Risk", 
     xlab="Month",lty=1,lwd=1,xlim=c(0,60), cex.lab=1.2,cex.main=1.2,main =c(paste("Risk of prostate cancer death"),paste("under elimination of competing events"))) #Death due to Prostate Cancer: Estimand (5)
lines(cutTimes, cumIncPlacebo3, type="s",col=2,ylim=c(0,1),lwd=1)
lines(cutTimes, cumIncTreatedIpw3, type="s",col=1,ylim=c(0,1),lwd=1,lty=2)
lines(cutTimes, cumIncPlaceboIpw3, type="s",col=2,ylim=c(0,1),lwd=1,lty=2)
legend("topleft", c("parametric g-formula (treatment)","parametric g-formula (placebo)","IPW (treatment)","IPW (placebo)"),
       col=c(1,2,1,2), lty=c(1,1,2,2),cex=1.2,pt.cex=1.2,lwd=2)
dev.off()

pdf("RiskCompetingEvent.pdf")
plot(cutTimes,cumIncTreated5comp, type="s",ylim=c(0,1), ylab="Risk", 
     xlab="Month",lty=1,lwd=1,xlim=c(0,60), cex.lab=1.2,cex.main=1.2, main = "Risk of other death (the competing event)") #Death due to other causes: Estimand (5)
lines(cutTimes, cumIncPlacebo5comp, type="s",col=2,ylim=c(0,1),lwd=1)
lines(cutTimes, cumIncTreatedIpw5_2comp, type="s",col=1,lwd=1,lty=2)
lines(cutTimes, cumIncPlaceboIpw5_2comp, type="s",col=2,lwd=1,lty=2)
lines(cutTimes, cumIncTreatedIpw5comp, type="s",col=1,ylim=c(0,1),lwd=1,lty=3)
lines(cutTimes, cumIncPlaceboIpw5comp, type="s",col=2,ylim=c(0,1),lwd=1,lty=3)
legend("topleft", c("parametric g-formula (treatment)","parametric g-formula (placebo)","IPW sub (treatment)","IPW sub (placebo)","IPW cs (treatment)","IPW cs (placebo)"),
       col=c(1,2,1,2,1,2), lty=c(1,1,2,2,3,3),cex=1.2,pt.cex=1.2,lwd=2)
dev.off()

pdf("RiskCompositeOutcome.pdf")
plot(cutTimes,cumIncTreated5+cumIncTreated5comp, type="s",ylim=c(0,1), ylab="Risk", 
     xlab="Month",lty=1,lwd=1,xlim=c(0,60), main = "Risk of death from any cause (composite outcome)") #Death due to other causes: Estimand (5)
lines(cutTimes, cumIncPlacebo5+cumIncPlacebo5comp, type="s",col=2,ylim=c(0,1),lwd=1)
lines(cutTimes, cumIncTreatedIpw5_2+cumIncTreatedIpw5_2comp, type="s",col=1,lwd=1,lty=2)
lines(cutTimes, cumIncPlaceboIpw5_2+cumIncPlaceboIpw5_2comp, type="s",col=2,lwd=1,lty=2)
lines(cutTimes, cumIncTreatedIpw5+cumIncTreatedIpw5comp, type="s",col=1,ylim=c(0,1),lwd=1,lty=3)
lines(cutTimes, cumIncPlaceboIpw5+cumIncPlaceboIpw5comp, type="s",col=2,ylim=c(0,1),lwd=1,lty=3)
legend("topleft", c("parametric g-formula (treatment)","parametric g-formula (placebo)","IPW sub (treatment)","IPW sub (placebo)","IPW cs (treatment)","IPW cs (placebo)"),
       col=c(1,2,1,2,1,2), lty=c(1,1,2,2,3,3),cex=1,pt.cex=1,lwd=2)
dev.off()


###### PLOTS for sanity check ######
