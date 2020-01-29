#R code for "A causal framework for classical statistical estimands in failure time settings with competing events"
#Manuscript authors: Jessica G. Young, Mats J. Stensrud, Eric J. Tchetgen Tchetgen, Miguel A. Hernan
#Code written by: Mats. J. Stensrud and Jessica G. Young

#############
###### This program implements an IPW estimator of risk of the competing event (here 'death from other causes') based on subdistribution hazards. 
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
prostRed$rx=abs(1-prostRed$temprx)#1 is DES, 0 is Placebo # JGY change so that DES is A=1 and placebo A=0
prostRed$eventType = as.integer(prostRed$eventType)-1 #0 is censoring, 1 is pdeath, 2 is odeath
prostRed$hgBinary <- prostRed$hg < 12
prostRed$ageCat <- cut2(prostRed$age,c(0,60,75,100))
prostRed$normalAct <- prostRed$pf == "normal activity"
prostRed$eventCens <- prostRed$eventType == 0
prostRed$eventProst <- prostRed$eventType == 1


prostRedIpwComp <- prostRed

#Make event times larger than follow-up for subjects experiencing the competing event. 
prostRedIpwComp$dtime_old <- prostRedIpwComp$dtime #Save the original event times
prostRedIpwComp$dtime <-  prostRedIpwComp$dtime * as.integer(prostRedIpwComp$eventType!=1)  + 100 * as.integer(prostRedIpwComp$eventType==1)


cutTimes <- c(0:59) # Follow-up for 5 years
prostRedIpwComp$Tstart = -0.01 #Needed for the long format

longProstRedIpwComp <- survSplit(data=prostRedIpwComp,cut=cutTimes,start="Tstart",end="dtime",event="allCause")
longProstRedIpwCompCens <- survSplit(data=prostRedIpwComp,cut=cutTimes,start="Tstart",end="dtime",event="eventCens")
longProstRedIpwComp$eventCens <- longProstRedIpwCompCens$eventCens
# Make column for prostate cancer mortality
longProstRedIpwComp$prostateDeath <- longProstRedIpwComp$allCause==1 & longProstRedIpwComp$eventType==1
longProstRedIpwComp$otherDeath <- longProstRedIpwComp$allCause==1 & longProstRedIpwComp$eventType==2
longProstRedIpwComp$otherDeath[longProstRedIpwComp$eventCens==1]<-NA
longProstRedIpwComp<-longProstRedIpwComp[longProstRedIpwComp$dtime<length(cutTimes),]


# Data processing: create powers of time to allow time dependent baseline hazard
longProstRedIpwComp$dtime2 <- longProstRedIpwComp$dtime * longProstRedIpwComp$dtime
longProstRedIpwComp$dtime3 <- longProstRedIpwComp$dtime2 * longProstRedIpwComp$dtime
longProstRedIpwComp$dtime4 <- longProstRedIpwComp$dtime3 * longProstRedIpwComp$dtime

# Data processing: create interactions with time and treatment
longProstRedIpwComp$rx1 <- longProstRedIpwComp$dtime * longProstRedIpwComp$rx
longProstRedIpwComp$rx2 <- longProstRedIpwComp$dtime2 * longProstRedIpwComp$rx
longProstRedIpwComp$rx3 <- longProstRedIpwComp$dtime3 * longProstRedIpwComp$rx
longProstRedIpwComp$rx4 <- longProstRedIpwComp$dtime4 * longProstRedIpwComp$rx

# Create 'baseline' - data collected at visit 0
baselineIpwComp <- longProstRedIpwComp[longProstRedIpwComp$dtime==0,]


# Find subjects (and event times) that do not contribute to the weigths
nonContributorsComp <- ( (longProstRedIpwComp$eventType==1) & (longProstRedIpwComp$dtime >= longProstRedIpwComp$dtime_old) ) | (longProstRedIpwComp$dtime<50)
#Fit censoring weights to the restricted sample not who did not experience competing event and month<=50
plrFitC_2comp <- glm(eventCens ~ normalAct+ageCat+hx 
                     + rx, data = longProstRedIpwComp[!nonContributorsComp,],family=binomial())

# Create weights for IpwComp estimation
predC_2comp <- 1-predict(plrFitC_2comp, newdata = longProstRedIpwComp, type = 'response')
predC_2comp[nonContributorsComp] <- 1 # For nonContributors, the weight contribution is 1. 
cumPredC_2comp <- unlist(aggregate(predC_2comp~longProstRedIpwComp$patno,FUN = cumprod)$predC,use.names = F) #unlist because aggregate creates list
IpwCompEstimand5_2comp <- 1/cumPredC_2comp


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
treatedIpw_hazardOcomp <- nonParametricCumHaz(IpwCompEstimand5_2comp, inputdata=longProstRedIpwComp, grp=1, outcomeProstate=FALSE)
treatedIpw_hazardPcomp <- rep(0,length.out=(length(cutTimes)))
cumIncTreatedIpw5_2comp <- nonParametricCumInc(treatedIpw_hazardOcomp,treatedIpw_hazardPcomp)


# Second for placebo
placeboIpw_hazardOcomp <- nonParametricCumHaz(IpwCompEstimand5_2comp, inputdata=longProstRedIpwComp, grp=0, outcomeProstate=FALSE)
placeboIpw_hazardPcomp <- rep(0,length.out=(length(cutTimes)))
cumIncPlaceboIpw5_2comp <- nonParametricCumInc(placeboIpw_hazardOcomp,placeboIpw_hazardPcomp)

IPWsubtotRRcomp<-cumIncTreatedIpw5_2comp[length(cutTimes)]/cumIncPlaceboIpw5_2comp[length(cutTimes)]
IPWsubtotRDcomp<-cumIncTreatedIpw5_2comp[length(cutTimes)]-cumIncPlaceboIpw5_2comp[length(cutTimes)]
print(IPWsubtotRRcomp)
print(IPWsubtotRDcomp)