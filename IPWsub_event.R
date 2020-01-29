#R code for "A causal framework for classical statistical estimands in failure time settings with competing events"
#Manuscript authors: Jessica G. Young, Mats J. Stensrud, Eric J. Tchetgen Tchetgen, Miguel A. Hernan
#Code written by: Mats. J. Stensrud and Jessica G. Young


#############
###### This program implements an estimator of the risk of the event of interest without elimination of competing events based on estimating subdistribution hazards (Here the outcome is prostate death). 
#############
setwd("your directory")
library("Hmisc")
prostate <- read.csv("prostate.csv")
str(prostate)
###########
###### Data processing
###########

# Create a binary variable to indicate all-cause death
prostate$allCause <- prostate$status != "alive"

# Create a variable with 3 levels indicating the cause of death 
prostate$eventType <- prostate$status
levels(prostate$eventType) <- list(alive="alive",pdeath="dead - prostatic ca",odeath = c(levels(prostate$eventType)[c(2:5,7:10)]))

# Data processing continues: Reduced data to only include high dose DES (A=1) and placebo (A=0), 
prostRed <- prostate[prostate$rx %in% levels(prostate$rx)[3:4],]
prostRed$temprx = as.integer(prostRed$rx)-3 
prostRed$rx=abs(1-prostRed$temprx)#1 is DES, 0 is Placebo 
prostRed$eventType = as.integer(prostRed$eventType)-1 #0 is censoring, 1 is pdeath, 2 is odeath
prostRed$hgBinary <- prostRed$hg < 12
prostRed$ageCat <- cut2(prostRed$age,c(0,60,75,100))
prostRed$normalAct <- prostRed$pf == "normal activity"
prostRed$eventCens <- prostRed$eventType == 0
prostRed$eventProst <- prostRed$eventType == 1



prostRedIpw <- prostRed

prostRedIpw$dtime_old <- prostRedIpw$dtime #Save the original event times
#Make event times larger than follow-up for subjects experiencing the competing event. 
prostRedIpw$dtime <-  prostRedIpw$dtime * as.integer(prostRedIpw$eventType!=2)  + 100 * as.integer(prostRedIpw$eventType==2)


cutTimes <- c(0:59) 
prostRedIpw$Tstart = -0.01 #Needed for the long format


longProstRedIpw <- survSplit(data=prostRedIpw,cut=cutTimes,start="Tstart",end="dtime",event="allCause")
longProstRedIpwCens <- survSplit(data=prostRedIpw,cut=cutTimes,start="Tstart",end="dtime",event="eventCens")
longProstRedIpw$eventCens <- longProstRedIpwCens$eventCens

# Make column for prostate cancer mortality
longProstRedIpw$prostateDeath <- longProstRedIpw$allCause==1 & longProstRedIpw$eventType==1
longProstRedIpw$otherDeath <- longProstRedIpw$allCause==1 & longProstRedIpw$eventType==2
longProstRedIpw$prostateDeath[longProstRedIpw$eventCens==1]<-NA
longProstRedIpw<-longProstRedIpw[longProstRedIpw$dtime<length(cutTimes),]
# Data processing: create powers of time to allow time dependent baseline hazard for censoring
longProstRedIpw$dtime2 <- longProstRedIpw$dtime * longProstRedIpw$dtime
longProstRedIpw$dtime3 <- longProstRedIpw$dtime2 * longProstRedIpw$dtime
longProstRedIpw$dtime4 <- longProstRedIpw$dtime3 * longProstRedIpw$dtime


# Data processing: create time interactions with treatment indicator
longProstRedIpw$rx1 <- longProstRedIpw$dtime * longProstRedIpw$rx
longProstRedIpw$rx2 <- longProstRedIpw$dtime2 * longProstRedIpw$rx
longProstRedIpw$rx3 <- longProstRedIpw$dtime3 * longProstRedIpw$rx
longProstRedIpw$rx4 <- longProstRedIpw$dtime4 * longProstRedIpw$rx

# Create 'baseline' - data collected at visit 0
baselineIpw <- longProstRedIpw[longProstRedIpw$dtime==0,]

# Find person-time records that do not contribute to the censoring weight models (probability of being censored for these is 0)
nonContributors <-  ((longProstRedIpw$eventType==2) & (longProstRedIpw$dtime >= longProstRedIpw$dtime_old)) | (longProstRedIpw$dtime<50)

#Fit censoring weights
plrFitC_2 <- glm(eventCens ~ normalAct+ageCat+hx 
                 + rx, data = longProstRedIpw[!nonContributors,],family=binomial())

## Create weights for Ipw estimation
predC_2 <- 1-predict(plrFitC_2, newdata = longProstRedIpw, type = 'response')
predC_2[nonContributors] <- 1 # For nonContributors, the weight contribution is 1. 
cumPredC_2 <- unlist(aggregate(predC_2~longProstRedIpw$patno,FUN = cumprod)$predC,use.names = F) #unlist because aggregate creates list
IpwEstimand5_2 <- 1/cumPredC_2
summary(IpwEstimand5_2)

nonParametricCumHaz <- function(weightVector, inputdata, grp, outcomeProstate=TRUE){
  outputHazards <- rep(NA, length.out=length(cutTimes))
  counter <- 1 
  for(i in cutTimes){
    if(outcomeProstate){
      indices <- inputdata$dtime==i & inputdata$rx == grp & inputdata$eventCens==0  #JGY & inputdata$otherDeath==0 I think need to take this out because want to keep in denominator those with competing event
      eventIndicator <- indices & inputdata$prostateDeath==1 
    }else{
      indices <- inputdata$dtime==i & inputdata$rx == grp & inputdata$eventCens==0
      eventIndicator <- indices & inputdata$otherDeath==1 
    }
    outputHazards[counter] <- sum(weightVector[eventIndicator]) / sum(weightVector[indices])
    counter <- counter+1
  }
  return(outputHazards)
}

nonParametricCumInc <- function(hazard1,hazard2,competing=FALSE){
  inc <- rep(NA, length.out=length(cutTimes))
  cumulativeSurvival <- c(1, cumprod( (1-hazard1) * (1-hazard2) ))
  counter <- 1 
  for(i in 1:length(cutTimes)){
    if(!competing){
      inc[i] <- hazard1[i] * (1-hazard2[i]) * cumulativeSurvival[i]
    }else{
      inc[i] <- hazard1[i] * cumulativeSurvival[i]
    }
  }
  cumInc <- cumsum(inc)
  return(cumInc)
}




#### Created dataset to predict hazards 

# Turn into conditional risk of event in each time interval
treatedIpw_hazardP <- nonParametricCumHaz(IpwEstimand5_2, inputdata=longProstRedIpw, grp=1, outcomeProstate=TRUE) 
treatedIpw_hazardO <- rep(0,length.out=(length(cutTimes)))
cumIncTreatedIpw5_2 <- nonParametricCumInc(treatedIpw_hazardP,treatedIpw_hazardO)


# Second for placebo
placeboIpw_hazardP <- nonParametricCumHaz(IpwEstimand5_2, inputdata=longProstRedIpw, grp=0, outcomeProstate=TRUE) 
placeboIpw_hazardO <- rep(0,length.out=(length(cutTimes)))
cumIncPlaceboIpw5_2 <- nonParametricCumInc(placeboIpw_hazardP,placeboIpw_hazardO)

IPWsubtotRR<-cumIncTreatedIpw5_2[length(cutTimes)]/cumIncPlaceboIpw5_2[length(cutTimes)] #JGY need ci's
IPWsubtotRD<-cumIncTreatedIpw5_2[length(cutTimes)]-cumIncPlaceboIpw5_2[length(cutTimes)] #JGY need ci's
print(IPWsubtotRR)
print(IPWsubtotRD)

