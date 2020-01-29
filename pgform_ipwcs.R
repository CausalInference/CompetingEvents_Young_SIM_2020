#R code for "A causal framework for classical statistical estimands in failure time settings with competing events"
#Manuscript authors: Jessica G. Young, Mats J. Stensrud, Eric J. Tchetgen Tchetgen, Miguel A. Hernan
#Code written by: Mats. J. Stensrud and Jessica G. Young

############ This program implements parametric g-formula estimators of total and controlled direct effects using pooled over time logistic models
############# as well as IPW estimators based on weighted estimates of the cause-specific hazards of the event of interest and competing events
############## (equivalently the "competing event conditioned hazard of the event of interest" and the "hazard of the competing event").

#The dataset is obtained from http://biostat.mc.vanderbilt.edu/DataSets".
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
# Convert prostRed to long format 
cutTimes <- c(0:59) 
prostRed$Tstart = -0.01 #Needed for the long format
longProstRed <- survSplit(data=prostRed,cut=cutTimes,start="Tstart",end="dtime",event="allCause")
longProstRedCens <- survSplit(data=prostRed,cut=cutTimes,start="Tstart",end="dtime",event="eventCens")
longProstRed$eventCens <- longProstRedCens$eventCens
# Make column for prostate cancer mortality
longProstRed$prostateDeath <- longProstRed$allCause==1 & longProstRed$eventType==1
longProstRed$otherDeath <- longProstRed$allCause==1 & longProstRed$eventType==2

#To be sure that right risk sets are used in models for this data structure without explicit conditioning, set "future values" of D and Y to missing on a given line if C=1 on that line and set future values of Y to missing on a line if D=1 on that line (respects "order" C,D,Y in each interval) -- nonissue when there are no "ties" but with month long intervals there might be.
longProstRed$prostateDeath[longProstRed$eventCens==1]<-NA
longProstRed$otherDeath[longProstRed$eventCens==1]<-NA
longProstRed$prostateDeath[longProstRed$otherDeath==1]<-NA

#Restrict input data set to records with dtime<K+1 for fitting pooled over time models
longProstRed<-longProstRed[longProstRed$dtime<length(cutTimes),]


# create powers of time to allow time dependent baseline hazard
longProstRed$dtime2 <- longProstRed$dtime * longProstRed$dtime
longProstRed$dtime3 <- longProstRed$dtime2 * longProstRed$dtime

# create time dependent coefficients for the hazard of prostate death
longProstRed$rx1 <- longProstRed$dtime * longProstRed$rx
longProstRed$rx2 <- longProstRed$dtime2 * longProstRed$rx
longProstRed$rx3 <- longProstRed$dtime3 * longProstRed$rx

# Create 'baseline' - data collected at visit 0
baseline <- longProstRed[longProstRed$dtime==0,]

# Number of subjects
n <- length(unique(longProstRed$patno))

############
###### Fit conditional pooled logistic regression models later used for parametric g-formula and/or ipw
####### Note all models are pooled over all patients
############
plrFitP <- glm(prostateDeath ~ dtime + dtime2 + dtime3 + normalAct+ageCat+hx+hgBinary 
               + rx + rx1 + rx2 , data = longProstRed,family=binomial())
plrFitO <- glm(otherDeath ~ dtime + dtime2  + normalAct+ageCat+hx+hgBinary 
                    + rx , data = longProstRed,family=binomial())


temp50<-longProstRed[longProstRed$dtime>50,]
plrFitC <-  plrFitC <- glm(eventCens ~  normalAct+ageCat+hx 
                                + rx , data = temp50,family=binomial())

# Fit models for the numerators in the stabilized weights MJS:
plrFitCnum <- glm(eventCens ~ rx , data = temp50,family=binomial())  
plrFitOnum <- glm(otherDeath ~ dtime + dtime2 + rx, data = longProstRed,family=binomial())

############
###### Create expanded datasets to obtain parametric g-formula estimates 
############

# Create expanded data so that hazards given A=1 can be computed for all subjects (even those with A=0) -- note in an RCT this is not necessary, could stratify models above and this step by treatment arm
# Expand baseline so it contains a visit at each time point for every individual
# where the baseline information has been carried forward at each time
treated <- baseline[rep(1:n,each=length(cutTimes)),] #One row for each interval k+1 for k=0,...,K
treated$dtime <- rep(cutTimes,n)
treated$rx <-1 
# Re-create the additional interaction terms. 
treated$dtime2 <- treated$dtime * treated$dtime
treated$dtime3 <- treated$dtime2 * treated$dtime
treated$dtime4 <- treated$dtime3 * treated$dtime

treated$rx1 <- treated$dtime * treated$rx
treated$rx2 <- treated$dtime2 * treated$rx
treated$rx3 <- treated$dtime3 * treated$rx
treated$rx4 <- treated$dtime4 * treated$rx

# Estimate conditional discrete hazards (cause specific for each cause of death/competing event conditioned for event and hazard of competing event) for each subject in each time interval
treated$hazardP <- predict(plrFitP, newdata = treated, type = 'response') 
treated$hazardO <- predict(plrFitO, newdata = treated, type = 'response')
treated$s <- (1-treated$hazardP) * (1-treated$hazardO)
# sum(treated$hazardO < 0) sanity check

#New data for risk under elimination of competing events ("3" and "5" suffixes are relic from equation numbering in early version of manuscript, "3" is for elimination and "5" without)
treated3 <- treated
treated3$hazardO <- 0
treated3$s <- (1-treated3$hazardP) * (1-treated3$hazardO)


# Make analogous dataset for placebo
placebo <- baseline[rep(1:n,each=length(cutTimes)),] #One row for each time
placebo$dtime <- rep(cutTimes,n)
placebo$rx <- 0 

# Re-create the additional interaction terms. 
placebo$dtime2 <- placebo$dtime * placebo$dtime
placebo$dtime3 <- placebo$dtime2 * placebo$dtime
placebo$dtime4 <- placebo$dtime3 * placebo$dtime

placebo$rx1 <- placebo$dtime * placebo$rx
placebo$rx2 <- placebo$dtime2 * placebo$rx
placebo$rx3 <- placebo$dtime3 * placebo$rx
placebo$rx4 <- placebo$dtime4 * placebo$rx

# Estimate conditional discrete hazard for each subject in each time interval
placebo$hazardP <- predict(plrFitP, newdata = placebo, type = 'response') 
placebo$hazardO <- predict(plrFitO, newdata = placebo, type = 'response')
placebo$s <- (1-placebo$hazardP) * (1-placebo$hazardO)

# new data for risk under elimination of competing events
placebo3 <- placebo
placebo3$hazardO <- 0
placebo3$s <- (1-placebo3$hazardP) * (1-placebo3$hazardO)

# Utility function for parametric g-formula estimators
calculateCumInc <- function(inputData,timepts=cutTimes,competing=FALSE){
  cumulativeIncidence <- matrix(NA,ncol = length(unique(inputData$patno)),nrow = length(timepts))
  # Insert the event probabilities at the first time interval: 
  # Needs to account for "temporal order" Dk+1 before Yk+1
  if(!competing) cumulativeIncidence[1,] <- inputData[inputData$dtime==0,]$hazardP*(1-inputData[inputData$dtime==0,]$hazardO)
  else cumulativeIncidence[1,] <- inputData[inputData$dtime==0,]$hazardO
  # Create a matrix compatible with 'cumulativeIncidence' with survival probabilities at each time 
  survivalProb <- t(aggregate(s~patno,data = inputData, FUN = cumprod)$s) #split the long data into subsets per subject
  for(i in 2:length(timepts)){
    subInputDataP <- inputData[inputData$dtime==(i-1),]$hazardP #OBS: dtime starts at 0
    subInputDataO <- inputData[inputData$dtime==(i-1),]$hazardO #OBS: dtime starts at 0 
    if(!competing) cumulativeIncidence[i,] <- subInputDataP * (1-subInputDataO) *  survivalProb[(i-1),] # OBS: survivalProb entry i is time point i-1
    else cumulativeIncidence[i,] <- subInputDataO * survivalProb[(i-1),]
  }
  meanCumulativeIncidence <- rowMeans(apply(cumulativeIncidence, MARGIN = 2,cumsum)) 
  return(meanCumulativeIncidence)
}



###########
###### Calculate risks and corresponding treatment effects (RR/RD) by end of follow-up using parametric g-formula 
###########
#parametric g-formula estimate of risk without elimination of competing events
cumIncTreated5 <- calculateCumInc(treated)#a=1
cumIncPlacebo5 <- calculateCumInc(placebo)#a=0
#parametric g-formula estimate of the total effect
gcomptotrr<-cumIncTreated5[length(cutTimes)]/cumIncPlacebo5[length(cutTimes)] #JGY need ci's on this, I think let's skip composite outcome this is all enough
gcomptotrd<-cumIncTreated5[length(cutTimes)]-cumIncPlacebo5[length(cutTimes)] #JGY need ci's on this
print(gcomptotrr)
print(gcomptotrd)


#parametric g-formula estimate of risk under elimination of competing events
cumIncTreated3 <- calculateCumInc(treated3)#a=1
cumIncPlacebo3 <- calculateCumInc(placebo3)#a=0
#gcomp estimate of the controlled direct effect
gcompcderr<-cumIncTreated3[length(cutTimes)]/cumIncPlacebo3[length(cutTimes)] #JGY need ci's on this
gcompcderd<-cumIncTreated3[length(cutTimes)]-cumIncPlacebo3[length(cutTimes)] #JGY need ci's on this
print(gcompcderr)
print(gcompcderd)

# parametric g-formula estimate of the risk of the competing event itself
cumIncTreated5comp <- calculateCumInc(treated, competing = TRUE)#a=1
cumIncPlacebo5comp <- calculateCumInc(placebo, competing = TRUE)#a=0
#(Total) effect on risk of competing event
gcompcomprr<-cumIncTreated5comp[length(cutTimes)]/cumIncPlacebo5comp[length(cutTimes)] #JGY need ci's on this
gcompcomprd<-cumIncTreated5comp[length(cutTimes)]-cumIncPlacebo5comp[length(cutTimes)] #JGY need ci's on this
print(gcompcomprr)
print(gcompcomprd)

##############
###### Find IPW estimates
##############

# Create unstabilized and stablized weights for IPW estimators
longProstRed$predC<-rep(1,dim(longProstRed)[1]) # start by setting all to 1
longProstRed$predC[longProstRed$dtime>50] <- 1-predict(plrFitC, newdata = temp50, type = 'response') #replace with model-based values for records with dtime>50
longProstRed$predCnum <-rep(1,dim(longProstRed)[1]) #JGY start by setting values on all rows to 1
longProstRed$predCnum [longProstRed$dtime>50] <- 1-predict(plrFitCnum, newdata = temp50, type = 'response') #replace with model-based values for records with dtime>50
predC<-longProstRed$predC 
predCnum <- longProstRed$predCnum 
predO <- 1-predict(plrFitO, newdata = longProstRed, type = 'response')
predOnum <- 1-predict(plrFitOnum, newdata = longProstRed, type = 'response') 

cumPredC <- unlist(aggregate(predC~longProstRed$patno,FUN = cumprod)$predC,use.names = F) #unlist because aggregate creates list
cumPredO <- unlist(aggregate(predO~longProstRed$patno,FUN = cumprod)$predO,use.names = F) #unlist because aggregate creates list
cumPredCnum <- unlist(aggregate(predCnum~longProstRed$patno,FUN = cumprod)$predC,use.names = F) # MJS
cumPredOnum <- unlist(aggregate(predOnum~longProstRed$patno,FUN = cumprod)$predO,use.names = F) # MJS

ipwEstimand3 <- 1/cumPredC * 1/cumPredO
ipwEstimand5 <- 1/cumPredC
ipwEstimand3stab <- cumPredCnum / cumPredC * cumPredOnum/cumPredO # MJS
ipwEstimand5stab <- cumPredCnum / cumPredC # MJS
summary(ipwEstimand3)
summary(ipwEstimand5)
summary(ipwEstimand3stab)
summary(ipwEstimand5stab)

summary(ipwEstimand3/ipwEstimand3stab)

#JGY: large weights for estimand3, it would be nice to try stabilized for this one if there is time.  would mean fitting a comparable model above for C and D depending only on dtime and rx

#utility functions for IPW estimators
nonParametricCumHaz <- function(weightVector, inputdata, grp, outcomeProstate=TRUE){
  outputHazards <- rep(NA, length.out=length(cutTimes))
  counter <- 1 
  for(i in cutTimes){
    if(outcomeProstate){
      indices <- inputdata$dtime==i & inputdata$rx == grp & inputdata$eventCens==0 & inputdata$otherDeath==0 
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




# IPW estimate of the risk under elimination of competing events and controlled direct effect on event of interest
TreatedIpw3_hazardP <- nonParametricCumHaz(ipwEstimand3, inputdata=longProstRed, grp=1, outcomeProstate=TRUE)
TreatedIpw3_hazardO <- rep(0,length.out=(length(cutTimes)))
cumIncTreatedIpw3 <- nonParametricCumInc(TreatedIpw3_hazardP,TreatedIpw3_hazardO)

PlaceboIpw3_hazardP <- nonParametricCumHaz(ipwEstimand3, inputdata=longProstRed, grp=0, outcomeProstate=TRUE)
PlaceboIpw3_hazardO <- rep(0,length.out=(length(cutTimes)))
cumIncPlaceboIpw3 <- nonParametricCumInc(PlaceboIpw3_hazardP,PlaceboIpw3_hazardO)

TreatedIpw3_hazardPs <- nonParametricCumHaz(ipwEstimand3stab, inputdata=longProstRed, grp=1, outcomeProstate=TRUE)
TreatedIpw3_hazardOs <- rep(0,length.out=(length(cutTimes)))
cumIncTreatedIpw3s <- nonParametricCumInc(TreatedIpw3_hazardPs,TreatedIpw3_hazardOs)

PlaceboIpw3_hazardPs <- nonParametricCumHaz(ipwEstimand3stab, inputdata=longProstRed, grp=0, outcomeProstate=TRUE)
PlaceboIpw3_hazardOs <- rep(0,length.out=(length(cutTimes)))
cumIncPlaceboIpw3s <- nonParametricCumInc(PlaceboIpw3_hazardPs,PlaceboIpw3_hazardOs)

ipwcderr<-cumIncTreatedIpw3[length(cutTimes)]/cumIncPlaceboIpw3[length(cutTimes)]#JGY need ci's for this
ipwcderd<-cumIncTreatedIpw3[length(cutTimes)]-cumIncPlaceboIpw3[length(cutTimes)]# JGY need ci's for this
print(ipwcderr)
print(ipwcderd)

# IPW  estimate of the risk without elimination of competing events and total effect on event of interest
treatedIpw5_hazardP <- nonParametricCumHaz(ipwEstimand5, inputdata=longProstRed, grp=1, outcomeProstate=TRUE)
treatedIpw5_hazardO <- nonParametricCumHaz(ipwEstimand5, inputdata=longProstRed, grp=1, outcomeProstate=FALSE)
cumIncTreatedIpw5 <- nonParametricCumInc(treatedIpw5_hazardP,treatedIpw5_hazardO)

placeboIpw5_hazardP <- nonParametricCumHaz(ipwEstimand5, inputdata=longProstRed, grp=0, outcomeProstate=TRUE)
placeboIpw5_hazardO <- nonParametricCumHaz(ipwEstimand5, inputdata=longProstRed, grp=0, outcomeProstate=FALSE)
cumIncPlaceboIpw5 <- nonParametricCumInc(placeboIpw5_hazardP, placeboIpw5_hazardO)

ipwtotrr<-cumIncTreatedIpw5[length(cutTimes)]/cumIncPlaceboIpw5[length(cutTimes)]#JGY need ci's for this
ipwtotrd<-cumIncTreatedIpw5[length(cutTimes)]-cumIncPlaceboIpw5[length(cutTimes)]# JGY need ci's for this
print(ipwtotrr)
print(ipwtotrd)

#IPW estimate of the risk of the competing event and (total) effect on competing event
cumIncTreatedIpw5comp <- nonParametricCumInc(treatedIpw5_hazardO,treatedIpw5_hazardP, competing = TRUE)
cumIncPlaceboIpw5comp <- nonParametricCumInc(placeboIpw5_hazardO, placeboIpw5_hazardP, competing = TRUE)

ipwcomprr<-cumIncTreatedIpw5comp[length(cutTimes)]/cumIncPlaceboIpw5comp[length(cutTimes)]#JGY need ci's for this
ipwcomprd<-cumIncTreatedIpw5comp[length(cutTimes)]-cumIncPlaceboIpw5comp[length(cutTimes)]# JGY need ci's for this
print(ipwcomprr)
print(ipwcomprd)

############
###### We have created separate programs for alternative IPW estimators of total effects based on estimates of the subdistribution hazard 
############
