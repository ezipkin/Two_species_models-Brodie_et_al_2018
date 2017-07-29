##Code to run the two species models
##Written seperately for each species pair
#The data for species 1 are in data0
#The data for species 2 are in data1

#Load the library to run the BUGS models
library(R2WinBUGS)

#--------------------------------- CIVET ------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))

data0<-read.csv('civet0_wide.csv', header=TRUE, sep=",")
data1<-read.csv('civet1_wide.csv', header=TRUE, sep=",")
camera<-read.csv('sampcovs_camhours_civet.csv', header=TRUE, sep=",")
studyday<-read.csv('sampcovs_studyday_civet.csv', header=TRUE, sep=",")
covars<-read.csv('sitecovs_civet.csv', header=TRUE, sep=",")
elev=covars$elev_avg200
slope=covars$elev_sd200

rownames(data0)=data0[,1]
data0=data0[,-1]
data0=as.matrix(data0)
rownames(data1)=data1[,1]
data1=data1[,-1]
data1=as.matrix(data1)

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
xx = 30 # used 64 in the original submission but Reviewer1 wanted a shorter sampling period (here 5 months)
data0 = data0[,1:xx]
data1 = data1[,1:xx]

# Info for histograms of detections
length(which(data0[]==1))
length(which(data0[]==2))
length(which(data0[]==3))
length(which(data0[]==4))
length(which(data0[]>=5))
length(which(data1[]==1))
length(which(data1[]==2))
length(which(data1[]==3))
length(which(data1[]==4))
length(which(data1[]>=5))

# Counts of detections for Table 1
num0 <- sum(data0, na.rm=T)
num1 <- sum(data1, na.rm=T)

# How many sightings for each species
average.site.count.0 = round(apply(data0, 1, mean, na.rm=TRUE), digits=2)
average.site.count.1 = round(apply(data1, 1, mean, na.rm=TRUE), digits=2)

################## SAMPLING COVARIATES #######################################
rownames(camera)=camera[,1]
camera=camera[,-1]
rownames(studyday)=studyday[,1]
studyday=studyday[,-1]

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
cam = camera[,1:xx]
cam = as.matrix(cam)
day = studyday[,1:xx]
day = as.matrix(day)

#Standardize day 
mday = mean(as.vector(day),na.rm=T)
sdday = sd(as.vector(day),na.rm=T)
day1 <- (day-mday)/sdday 
for (i in 1:dim(day1)[1]) {
  a <- which(day[i,]>0)
  day1[i,-a]=0}

#Standardize camera hours 
mcam=mean(as.vector(cam),na.rm=T)
sdcam=sd(as.vector(cam),na.rm=T)
cam1<-(cam-mcam)/sdcam
for (i in 1:dim(cam1)[1]) {
  a <- which(cam[i,]>0)
  cam1[i,-a]=0}

########################## ABUNDANCE COVARIATES ########################
# Check that the covariate data are in the same order as the abundance data
#a=pmatch(dimnames(data0)[[1]],covars$station)
#a-nrow(covars)

# Elevation & slope standardized
melev=mean(elev,na.rm=T)
sdelev=sd(elev,na.rm=T)
elev1<-(elev-melev)/sdelev 
mslope=mean(slope,na.rm=T)
sdslope=sd(slope,na.rm=T)
slope1<-(slope-mslope)/sdslope

# Elevation & slope not standardized
#elev1<-elev
#slope1<-slope

#Water
river=covars$water
mriver=mean(river,na.rm=T)
sdriver=sd(river,na.rm=T)
river1<-(river-mriver)/sdriver 

#Logged
log=covars$logged
#log=covars$newlog

#trail
trail=covars$trail

#Number of sites
nsite=dim(data0)[1]

#Pull out area for random effect
area1 = (as.character(covars$area))
uarea = sort(unique(area1))
narea = length(uarea)

area=rep(NA,length(area1))
for (i in 1:length(area)) {
  area[i] = pmatch(area1[i],uarea)
}

########################### MODEL #######################################
sink("covar3.txt")
cat("
model{
  # Priors
  for (i in 1:2) {
  a0[i] ~ dnorm(0,0.01)
  a1[i] ~ dnorm(0,0.01)
  a11[i] ~ dnorm(0,0.01)
  a2[i] ~ dnorm(0,0.01)
  a3[i] ~ dnorm(0,0.01)
  a4[i] ~ dnorm(0,0.01)
  a5[i] ~ dnorm(0,0.01)
  b0[i] ~ dnorm(0,0.01)
  b1[i] ~ dnorm(0,0.01)
  b2[i] ~ dnorm(0,0.01)
  b3[i] ~ dnorm(0,0.01)
  b4[i] ~ dnorm(0,0.01)
}
  a6 ~ dnorm(0,0.01)

#var.a7 ~ dgamma(0.01,0.01)
#var.a8 ~ dgamma(0.01,0.01)

#var.a7 ~ dunif(0.01,5)
#var.a8 ~ dunif(0.01,10)

sigma.a7 ~ dunif(0,5)
var.a7 <- 1/(sigma.a7*sigma.a7)

sigma.a8 ~ dunif(0,10)
var.a8 <- 1/(sigma.a8*sigma.a8)

for (k in 1:narea) {
  a7[k] ~ dnorm(0,var.a7)
  a8[k] ~ dnorm(0,var.a8)
}

  # Likelihood
  # Ecological model for true abundance
  for (j in 1:nsite) {

    N0[j] ~ dpois(lambda0[j])

    log(lambda0[j]) <- a0[1] + a1[1]*slope1[j] + a11[1]*slope1[j]*slope1[j] + a2[1]*elev1[j] + a3[1]*elev1[j]*elev1[j] +
                      a4[1]*river1[j] + a5[1]*log[j] + a6*N1[j] + a7[area[j]]
    N1[j] ~ dpois(lambda1[j])
    log(lambda1[j]) <- a0[2] + a1[2]*slope1[j] + a11[2]*slope1[j]*slope1[j] + a2[2]*elev1[j] + a3[2]*elev1[j]*elev1[j] +
                      a4[2]*river1[j] + a5[2]*log[j] + a8[area[j]]
    # Observation model for replicated counts
    for (k in 1:nrep) {
      y0[j,k] ~ dbin(p0[j,k], N0[j])
        p0[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
      lp[j,k] <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
      
      y1[j,k] ~ dbin(p1[j,k], N1[j])
        p1[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
      llp[j,k] <- b0[2] + b1[2]*trail[j] + b2[2]*cam1[j,k] + b3[2]*day1[j,k] + b4[2]*log[j]

    #probit(p0[j,k]) <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
    } #k
  } #j

  #Derived parameters
  diff.slope <-  a1[1] - a1[2]
  diff.slope2 <-  a11[1] - a11[2]
  diff.elev <- a2[1] - a2[2]
  diff.elev2 <-  a3[1] - a3[2]
  diff.river <- a4[1] - a4[2]
  diff.log <-  a5[1] - a5[2]

} 
    ",fill = TRUE)
sink()

################################# BAYESIAN MODELING #############################
# Create the necessary arguments to run the bugs() command 
sp.data = list(y0=data0, y1=data1, nsite=nsite, nrep=dim(data0)[2], narea=narea, area=area,
               cam1=cam1, day1=day1, trail=trail, elev1=elev1, slope1=slope1, river1=river1, log=log) 

# Specify the parameters to be monitored
sp.params = list('a0', 'a1', 'a11', 'a2', 'a3', 'a4', 'a5', 'a6', 'b0', 'b1', 'b2', 'b3', 'b4', 'N0', 'N1', 'lambda0', 'lambda1',
	'diff.slope', 'diff.slope2', 'diff.elev', 'diff.elev2', 'diff.river', 'diff.log',
	'var.a7', 'var.a8', 'a7', 'a8')

# Specify the initial values
sp.inits = function() {
    list(a0=rnorm(2), a1=rnorm(2), a11=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2), a6=rnorm(1),   
         b0=rnorm(2), b1=rnorm(2), b2=rnorm(2), b3=rnorm(2), b4=rnorm(2), a7=rnorm(narea), a8=rnorm(narea),
         sigma.a7=runif(1,0,3), sigma.a8=runif(1,0,3),
         N0 = as.vector(apply(data0,1,max, na.rm=T) + 5),   
         N1 = as.vector(apply(data1,1,max, na.rm=T) + 5)
    )           
}

# MCMC settings (choose one)
ni <- 300; nt <- 2; nb <- 100; nc <- 3 # quick test
ni <- 20000; nt <- 10; nb <- 10000; nc <- 3 # examining parameter values
ni <- 100000; nt <- 20; nb <- 50000; nc <- 3 # publication-quality run

# CIVET
Sys.time()
out <- bugs(sp.data, sp.inits, sp.params, "covar3.txt", bugs.directory="c:/MyPrograms/WinBUGS14/",
	n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F)
save.image(file="workspace_civet.RData")
Sys.time()

#plot(out)
#out$summary  #Show mean, sd, CIs, second to last column is the Rhat, last column is the effect size



#--------------------------------- MACAQUE ------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))

data0<-read.csv('macaque0_wide.csv', header=TRUE, sep=",")
data1<-read.csv('macaque1_wide.csv', header=TRUE, sep=",")

# analyzing macaques by group rather than by individual because in large groups detections are non-indep
d0_first<-data0[,1]
d0_others<-data0[,-1]
d0_others[d0_others > 1] <- 1
dat0 <- cbind(d0_first,d0_others)
names(dat0)[1] <- "station"
data0 <- dat0
d1_first<-data1[,1]
d1_others<-data1[,-1]
d1_others[d1_others > 1] <- 1
dat1 <- cbind(d1_first,d1_others)
names(dat1)[1] <- "station"
data1 <- dat1

camera<-read.csv('sampcovs_camhours_maca.csv', header=TRUE, sep=",")
studyday<-read.csv('sampcovs_studyday_maca.csv', header=TRUE, sep=",")
covars<-read.csv('sitecovs_maca.csv', header=TRUE, sep=",")
elev=covars$elev_avg100
slope=covars$elev_sd100

rownames(data0)=data0[,1]
data0=data0[,-1]
data0=as.matrix(data0)
rownames(data1)=data1[,1]
data1=data1[,-1]
data1=as.matrix(data1)

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
xx = 30 # used 64 in the original submission but Reviewer1 wanted a shorter sampling period (here 5 months)
data0 = data0[,1:xx]
data1 = data1[,1:xx]

# Info for histograms of detections
length(which(data0[]==1))
length(which(data0[]==2))
length(which(data0[]==3))
length(which(data0[]==4))
length(which(data0[]>=5))
length(which(data1[]==1))
length(which(data1[]==2))
length(which(data1[]==3))
length(which(data1[]==4))
length(which(data1[]>=5))

# Counts of detections for Table 1
num0 <- sum(data0, na.rm=T)
num1 <- sum(data1, na.rm=T)

# How many sightings for each species
average.site.count.0 = round(apply(data0, 1, mean, na.rm=TRUE), digits=2)
average.site.count.1 = round(apply(data1, 1, mean, na.rm=TRUE), digits=2)


################## SAMPLING COVARIATES #######################################
rownames(camera)=camera[,1]
camera=camera[,-1]
rownames(studyday)=studyday[,1]
studyday=studyday[,-1]

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
cam = camera[,1:xx]
cam = as.matrix(cam)
day = studyday[,1:xx]
day = as.matrix(day)

#Standardize day 
mday = mean(as.vector(day),na.rm=T)
sdday = sd(as.vector(day),na.rm=T)
day1 <- (day-mday)/sdday 
for (i in 1:dim(day1)[1]) {
  a <- which(day[i,]>0)
  day1[i,-a]=0}

#Standardize camera hours 
mcam=mean(as.vector(cam),na.rm=T)
sdcam=sd(as.vector(cam),na.rm=T)
cam1<-(cam-mcam)/sdcam
for (i in 1:dim(cam1)[1]) {
  a <- which(cam[i,]>0)
  cam1[i,-a]=0}

########################## ABUNDANCE COVARIATES ########################
# Check that the covariate data are in the same order as the abundance data
#a=pmatch(dimnames(data0)[[1]],covars$station)
#a-nrow(covars)

# Elevation & slope standardized
melev=mean(elev,na.rm=T)
sdelev=sd(elev,na.rm=T)
elev1<-(elev-melev)/sdelev 
mslope=mean(slope,na.rm=T)
sdslope=sd(slope,na.rm=T)
slope1<-(slope-mslope)/sdslope

# Elevation & slope not standardized
#elev1<-elev
#slope1<-slope

#Water
river=covars$water
mriver=mean(river,na.rm=T)
sdriver=sd(river,na.rm=T)
river1<-(river-mriver)/sdriver 

#Logged
log=covars$logged
#log=covars$newlog

#trail
trail=covars$trail

#Number of sites
nsite=dim(data0)[1]

#Pull out area for random effect
area1 = (as.character(covars$area))
uarea = sort(unique(area1))
narea = length(uarea)

area=rep(NA,length(area1))
for (i in 1:length(area)) {
  area[i] = pmatch(area1[i],uarea)
}

########################### MODEL #######################################
sink("covar3.txt")
cat("
model{
  # Priors
  for (i in 1:2) {
  a0[i] ~ dnorm(0,0.01)
  a1[i] ~ dnorm(0,0.01)
  a11[i] ~ dnorm(0,0.01)
  a2[i] ~ dnorm(0,0.01)
  a3[i] ~ dnorm(0,0.01)
  a4[i] ~ dnorm(0,0.01)
  a5[i] ~ dnorm(0,0.01)
  b0[i] ~ dnorm(0,0.01)
  b1[i] ~ dnorm(0,0.01)
  b2[i] ~ dnorm(0,0.01)
  b3[i] ~ dnorm(0,0.01)
  b4[i] ~ dnorm(0,0.01)
}
  a6 ~ dnorm(0,0.01)

#var.a7 ~ dgamma(0.01,0.01)
#var.a8 ~ dgamma(0.01,0.01)

#var.a7 ~ dunif(0.01,5)
#var.a8 ~ dunif(0.01,10)

sigma.a7 ~ dunif(0,5)
var.a7 <- 1/(sigma.a7*sigma.a7)

sigma.a8 ~ dunif(0,10)
var.a8 <- 1/(sigma.a8*sigma.a8)

for (k in 1:narea) {
  a7[k] ~ dnorm(0,var.a7)
  a8[k] ~ dnorm(0,var.a8)
}

  # Likelihood
  # Ecological model for true abundance
  for (j in 1:nsite) {

    N0[j] ~ dpois(lambda0[j])

    log(lambda0[j]) <- a0[1] + a1[1]*slope1[j] + a11[1]*slope1[j]*slope1[j] + a2[1]*elev1[j] + a3[1]*elev1[j]*elev1[j] +
                      a4[1]*river1[j] + a5[1]*log[j] + a6*N1[j] + a7[area[j]]
    N1[j] ~ dpois(lambda1[j])
    log(lambda1[j]) <- a0[2] + a1[2]*slope1[j] + a11[2]*slope1[j]*slope1[j] + a2[2]*elev1[j] + a3[2]*elev1[j]*elev1[j] +
                      a4[2]*river1[j] + a5[2]*log[j] + a8[area[j]]
    # Observation model for replicated counts
    for (k in 1:nrep) {
      y0[j,k] ~ dbin(p0[j,k], N0[j])
        p0[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
      lp[j,k] <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
      
      y1[j,k] ~ dbin(p1[j,k], N1[j])
        p1[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
      llp[j,k] <- b0[2] + b1[2]*trail[j] + b2[2]*cam1[j,k] + b3[2]*day1[j,k] + b4[2]*log[j]

    #probit(p0[j,k]) <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
    } #k
  } #j

  #Derived parameters
  diff.slope <-  a1[1] - a1[2]
  diff.slope2 <-  a11[1] - a11[2]
  diff.elev <- a2[1] - a2[2]
  diff.elev2 <-  a3[1] - a3[2]
  diff.river <- a4[1] - a4[2]
  diff.log <-  a5[1] - a5[2]

} 
    ",fill = TRUE)
sink()

################################# BAYESIAN MODELING #############################
# Create the necessary arguments to run the bugs() command 
sp.data = list(y0=data0, y1=data1, nsite=nsite, nrep=dim(data0)[2], narea=narea, area=area,
               cam1=cam1, day1=day1, trail=trail, elev1=elev1, slope1=slope1, river1=river1, log=log) 

# Specify the parameters to be monitored
sp.params = list('a0', 'a1', 'a11', 'a2', 'a3', 'a4', 'a5', 'a6', 'b0', 'b1', 'b2', 'b3', 'b4', 'N0', 'N1', 'lambda0', 'lambda1',
	'diff.slope', 'diff.slope2', 'diff.elev', 'diff.elev2', 'diff.river', 'diff.log',
	'var.a7', 'var.a8', 'a7', 'a8')

# Specify the initial values
sp.inits = function() {
    list(a0=rnorm(2), a1=rnorm(2), a11=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2), a6=rnorm(1),   
         b0=rnorm(2), b1=rnorm(2), b2=rnorm(2), b3=rnorm(2), b4=rnorm(2), a7=rnorm(narea), a8=rnorm(narea),
         sigma.a7=runif(1,0,3), sigma.a8=runif(1,0,3),
         N0 = as.vector(apply(data0,1,max, na.rm=T) + 5),   
         N1 = as.vector(apply(data1,1,max, na.rm=T) + 5)
    )           
}

# MCMC settings (choose one)
ni <- 300; nt <- 2; nb <- 100; nc <- 3 # quick test
ni <- 20000; nt <- 10; nb <- 10000; nc <- 3 # examining parameter values
ni <- 100000; nt <- 20; nb <- 50000; nc <- 3 # publication-quality run

# MACAQUE
Sys.time()
out <- bugs(sp.data, sp.inits, sp.params, "covar3.txt", bugs.directory="c:/MyPrograms/WinBUGS14/",
	n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F)
save.image(file="workspace_macaque.RData")
Sys.time()

#plot(out)
#out$summary  #Show mean, sd, CIs, second to last column is the Rhat, last column is the effect size



#--------------------------------- MUNTJAC ------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))

data0<-read.csv('muntjac0_wide.csv', header=TRUE, sep=",")
data1<-read.csv('muntjac1_wide.csv', header=TRUE, sep=",")
camera<-read.csv('sampcovs_camhours_munt.csv', header=TRUE, sep=",")
studyday<-read.csv('sampcovs_studyday_munt.csv', header=TRUE, sep=",")
covars<-read.csv('sitecovs_munt.csv', header=TRUE, sep=",")
elev=covars$elev_avg100
slope=covars$elev_sd100

rownames(data0)=data0[,1]
data0=data0[,-1]
data0=as.matrix(data0)
rownames(data1)=data1[,1]
data1=data1[,-1]
data1=as.matrix(data1)

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
xx = 30 # used 64 in the original submission but Reviewer1 wanted a shorter sampling period (here 5 months)
data0 = data0[,1:xx]
data1 = data1[,1:xx]

# Counts of detections for Table 1
num0 <- sum(data0, na.rm=T)
num1 <- sum(data1, na.rm=T)

# Info for histograms of detections
length(which(data0[]==1))
length(which(data0[]==2))
length(which(data0[]==3))
length(which(data0[]==4))
length(which(data0[]>=5))
length(which(data1[]==1))
length(which(data1[]==2))
length(which(data1[]==3))
length(which(data1[]==4))
length(which(data1[]>=5))

# How many sightings for each species
average.site.count.0 = round(apply(data0, 1, mean, na.rm=TRUE), digits=2)
average.site.count.1 = round(apply(data1, 1, mean, na.rm=TRUE), digits=2)


################## SAMPLING COVARIATES #######################
rownames(camera)=camera[,1]
camera=camera[,-1]
rownames(studyday)=studyday[,1]
studyday=studyday[,-1]

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
cam = camera[,1:xx]
cam = as.matrix(cam)
day = studyday[,1:xx]
day = as.matrix(day)

#Standardize day 
mday = mean(as.vector(day),na.rm=T)
sdday = sd(as.vector(day),na.rm=T)
day1 <- (day-mday)/sdday 
for (i in 1:dim(day1)[1]) {
  a <- which(day[i,]>0)
  day1[i,-a]=0}

#Standardize camera hours 
mcam=mean(as.vector(cam),na.rm=T)
sdcam=sd(as.vector(cam),na.rm=T)
cam1<-(cam-mcam)/sdcam
for (i in 1:dim(cam1)[1]) {
  a <- which(cam[i,]>0)
  cam1[i,-a]=0}

########################## ABUNDANCE COVARIATES #################
# Check that the covariate data are in the same order as the abundance data
#a=pmatch(dimnames(data0)[[1]],covars$station)
#a-nrow(covars)

# Elevation & slope standardized
melev=mean(elev,na.rm=T)
sdelev=sd(elev,na.rm=T)
elev1<-(elev-melev)/sdelev 
mslope=mean(slope,na.rm=T)
sdslope=sd(slope,na.rm=T)
slope1<-(slope-mslope)/sdslope

# Elevation & slope not standardized
#elev1<-elev
#slope1<-slope

#Water
river=covars$water
mriver=mean(river,na.rm=T)
sdriver=sd(river,na.rm=T)
river1<-(river-mriver)/sdriver 

#Logged
log=covars$logged
#log=covars$newlog

#trail
trail=covars$trail

#Number of sites
nsite=dim(data0)[1]

#Pull out area for random effect
area1 = (as.character(covars$area))
uarea = sort(unique(area1))
narea = length(uarea)

area=rep(NA,length(area1))
for (i in 1:length(area)) {
  area[i] = pmatch(area1[i],uarea)
}

########################### MODEL #######################################
sink("covar3.txt")
cat("
model{
  # Priors
  for (i in 1:2) {
  a0[i] ~ dnorm(0,0.01)
  a1[i] ~ dnorm(0,0.01)
  a11[i] ~ dnorm(0,0.01)
  a2[i] ~ dnorm(0,0.01)
  a3[i] ~ dnorm(0,0.01)
  a4[i] ~ dnorm(0,0.01)
  a5[i] ~ dnorm(0,0.01)
  b0[i] ~ dnorm(0,0.01)
  b1[i] ~ dnorm(0,0.01)
  b2[i] ~ dnorm(0,0.01)
  b3[i] ~ dnorm(0,0.01)
  b4[i] ~ dnorm(0,0.01)
}
  a6 ~ dnorm(0,0.01)

#var.a7 ~ dgamma(0.01,0.01)
#var.a8 ~ dgamma(0.01,0.01)

#var.a7 ~ dunif(0.01,5)
#var.a8 ~ dunif(0.01,10)

sigma.a7 ~ dunif(0,5)
var.a7 <- 1/(sigma.a7*sigma.a7)

sigma.a8 ~ dunif(0,10)
var.a8 <- 1/(sigma.a8*sigma.a8)

for (k in 1:narea) {
  a7[k] ~ dnorm(0,var.a7)
  a8[k] ~ dnorm(0,var.a8)
}

  # Likelihood
  # Ecological model for true abundance
  for (j in 1:nsite) {

    N0[j] ~ dpois(lambda0[j])

    log(lambda0[j]) <- a0[1] + a1[1]*slope1[j] + a11[1]*slope1[j]*slope1[j] + a2[1]*elev1[j] + a3[1]*elev1[j]*elev1[j] +
                      a4[1]*river1[j] + a5[1]*log[j] + a6*N1[j] + a7[area[j]]
    N1[j] ~ dpois(lambda1[j])
    log(lambda1[j]) <- a0[2] + a1[2]*slope1[j] + a11[2]*slope1[j]*slope1[j] + a2[2]*elev1[j] + a3[2]*elev1[j]*elev1[j] +
                      a4[2]*river1[j] + a5[2]*log[j] + a8[area[j]]
    # Observation model for replicated counts
    for (k in 1:nrep) {
      y0[j,k] ~ dbin(p0[j,k], N0[j])
        p0[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
      lp[j,k] <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
      
      y1[j,k] ~ dbin(p1[j,k], N1[j])
        p1[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
      llp[j,k] <- b0[2] + b1[2]*trail[j] + b2[2]*cam1[j,k] + b3[2]*day1[j,k] + b4[2]*log[j]

    #probit(p0[j,k]) <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
    } #k
  } #j

  #Derived parameters
  diff.slope <-  a1[1] - a1[2]
  diff.slope2 <-  a11[1] - a11[2]
  diff.elev <- a2[1] - a2[2]
  diff.elev2 <-  a3[1] - a3[2]
  diff.river <- a4[1] - a4[2]
  diff.log <-  a5[1] - a5[2]

} 
    ",fill = TRUE)
sink()

################################# BAYESIAN MODELING #############################
# Create the necessary arguments to run the bugs() command 
sp.data = list(y0=data0, y1=data1, nsite=nsite, nrep=dim(data0)[2], narea=narea, area=area,
               cam1=cam1, day1=day1, trail=trail, elev1=elev1, slope1=slope1, river1=river1, log=log) 

# Specify the parameters to be monitored
sp.params = list('a0', 'a1', 'a11', 'a2', 'a3', 'a4', 'a5', 'a6', 'b0', 'b1', 'b2', 'b3', 'b4', 'N0', 'N1', 'lambda0', 'lambda1',
	'diff.slope', 'diff.slope2', 'diff.elev', 'diff.elev2', 'diff.river', 'diff.log',
	'var.a7', 'var.a8', 'a7', 'a8')

# Specify the initial values
sp.inits = function() {
    list(a0=rnorm(2), a1=rnorm(2), a11=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2), a6=rnorm(1),   
         b0=rnorm(2), b1=rnorm(2), b2=rnorm(2), b3=rnorm(2), b4=rnorm(2), a7=rnorm(narea), a8=rnorm(narea),
         sigma.a7=runif(1,0,3), sigma.a8=runif(1,0,3),
         N0 = as.vector(apply(data0,1,max, na.rm=T) + 5),   
         N1 = as.vector(apply(data1,1,max, na.rm=T) + 5)
    )           
}

# MCMC settings (choose one)
ni <- 300; nt <- 2; nb <- 100; nc <- 3 # quick test
ni <- 20000; nt <- 10; nb <- 10000; nc <- 3 # examining parameter values
ni <- 100000; nt <- 20; nb <- 50000; nc <- 3 # publication-quality run

# MUNTJAC
Sys.time()
out <- bugs(sp.data, sp.inits, sp.params, "covar3.txt", bugs.directory="c:/MyPrograms/WinBUGS14/",
	n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F)
save.image(file="workspace_munt.RData")
Sys.time()

#plot(out)
#out$summary  #Show mean, sd, CIs, second to last column is the Rhat, last column is the effect size



#--------------------------------- PORCUPINE ------------------------------------------------------------------------------------
rm(list=ls(all=TRUE))

data0<-read.csv('porcupine0_wide.csv', header=TRUE, sep=",")
data1<-read.csv('porcupine1_wide.csv', header=TRUE, sep=",")
camera<-read.csv('sampcovs_camhours_porc.csv', header=TRUE, sep=",")
studyday<-read.csv('sampcovs_studyday_porc.csv', header=TRUE, sep=",")
covars<-read.csv('sitecovs_porc.csv', header=TRUE, sep=",")
elev=covars$elev_avg100
slope=covars$elev_sd100

# New to this iteration -analyzing porcupines by group rather than by individual
d0_first<-data0[,1]
d0_others<-data0[,-1]
d0_others[d0_others > 1] <- 1
dat0 <- cbind(d0_first,d0_others)
names(dat0)[1] <- "station"
data0 <- dat0
d1_first<-data1[,1]
d1_others<-data1[,-1]
d1_others[d1_others > 1] <- 1
dat1 <- cbind(d1_first,d1_others)
names(dat1)[1] <- "station"
data1 <- dat1

rownames(data0)=data0[,1]
data0=data0[,-1]
data0=as.matrix(data0)
rownames(data1)=data1[,1]
data1=data1[,-1]
data1=as.matrix(data1)

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
xx = 30 # used 64 in the original submission but Reviewer1 wanted a shorter sampling period (here 5 months)
data0 = data0[,1:xx]
data1 = data1[,1:xx]

# Counts of detections for Table 1
num0 <- sum(data0, na.rm=T)
num1 <- sum(data1, na.rm=T)

# Info for histograms of detections
length(which(data0[]==1))
length(which(data0[]==2))
length(which(data0[]==3))
length(which(data0[]==4))
length(which(data0[]>=5))
length(which(data1[]==1))
length(which(data1[]==2))
length(which(data1[]==3))
length(which(data1[]==4))
length(which(data1[]>=5))

# How many sightings for each species
average.site.count.0 = round(apply(data0, 1, mean, na.rm=TRUE), digits=2)
average.site.count.1 = round(apply(data1, 1, mean, na.rm=TRUE), digits=2)
 

################## SAMPLING COVARIATES ########################
rownames(camera)=camera[,1]
camera=camera[,-1]
rownames(studyday)=studyday[,1]
studyday=studyday[,-1]

#There are a lot of NAs in this data and detection appears very low
#For the time being I am going to truncate the data to the first xx reps
#It probably makes more sense to combine these reps
cam = camera[,1:xx]
cam = as.matrix(cam)
day = studyday[,1:xx]
day = as.matrix(day)

#Standardize day 
mday = mean(as.vector(day),na.rm=T)
sdday = sd(as.vector(day),na.rm=T)
day1 <- (day-mday)/sdday 
for (i in 1:dim(day1)[1]) {
  a <- which(day[i,]>0)
  day1[i,-a]=0}

#Standardize camera hours 
mcam=mean(as.vector(cam),na.rm=T)
sdcam=sd(as.vector(cam),na.rm=T)
cam1<-(cam-mcam)/sdcam
for (i in 1:dim(cam1)[1]) {
  a <- which(cam[i,]>0)
  cam1[i,-a]=0}

########################## ABUNDANCE COVARIATES ###############
# Check that the covariate data are in the same order as the abundance data
#a=pmatch(dimnames(data0)[[1]],covars$station)
#a-nrow(covars)

# Elevation & slope standardized
melev=mean(elev,na.rm=T)
sdelev=sd(elev,na.rm=T)
elev1<-(elev-melev)/sdelev 
mslope=mean(slope,na.rm=T)
sdslope=sd(slope,na.rm=T)
slope1<-(slope-mslope)/sdslope

# Elevation & slope not standardized
#elev1<-elev
#slope1<-slope

#Water
river=covars$water
mriver=mean(river,na.rm=T)
sdriver=sd(river,na.rm=T)
river1<-(river-mriver)/sdriver 

#Logged
log=covars$logged
#log=covars$newlog

#trail
trail=covars$trail

#Number of sites
nsite=dim(data0)[1]

#Pull out area for random effect
area1 = (as.character(covars$area))
uarea = sort(unique(area1))
narea = length(uarea)

area=rep(NA,length(area1))
for (i in 1:length(area)) {
  area[i] = pmatch(area1[i],uarea)
}

########################### MODEL #######################################
sink("covar3.txt")
cat("
model{
  # Priors
  for (i in 1:2) {
  a0[i] ~ dnorm(0,0.01)
  a1[i] ~ dnorm(0,0.01)
  a11[i] ~ dnorm(0,0.01)
  a2[i] ~ dnorm(0,0.01)
  a3[i] ~ dnorm(0,0.01)
  a4[i] ~ dnorm(0,0.01)
  a5[i] ~ dnorm(0,0.01)
  b0[i] ~ dnorm(0,0.01)
  b1[i] ~ dnorm(0,0.01)
  b2[i] ~ dnorm(0,0.01)
  b3[i] ~ dnorm(0,0.01)
  b4[i] ~ dnorm(0,0.01)
}
  a6 ~ dnorm(0,0.01)

#var.a7 ~ dgamma(0.01,0.01)
#var.a8 ~ dgamma(0.01,0.01)

#var.a7 ~ dunif(0.01,5)
#var.a8 ~ dunif(0.01,10)

sigma.a7 ~ dunif(0,5)
var.a7 <- 1/(sigma.a7*sigma.a7)

sigma.a8 ~ dunif(0,10)
var.a8 <- 1/(sigma.a8*sigma.a8)

for (k in 1:narea) {
  a7[k] ~ dnorm(0,var.a7)
  a8[k] ~ dnorm(0,var.a8)
}

  # Likelihood
  # Ecological model for true abundance
  for (j in 1:nsite) {

    N0[j] ~ dpois(lambda0[j])

    log(lambda0[j]) <- a0[1] + a1[1]*slope1[j] + a11[1]*slope1[j]*slope1[j] + a2[1]*elev1[j] + a3[1]*elev1[j]*elev1[j] +
                      a4[1]*river1[j] + a5[1]*log[j] + a6*N1[j] + a7[area[j]]
    N1[j] ~ dpois(lambda1[j])
    log(lambda1[j]) <- a0[2] + a1[2]*slope1[j] + a11[2]*slope1[j]*slope1[j] + a2[2]*elev1[j] + a3[2]*elev1[j]*elev1[j] +
                      a4[2]*river1[j] + a5[2]*log[j] + a8[area[j]]
    # Observation model for replicated counts
    for (k in 1:nrep) {
      y0[j,k] ~ dbin(p0[j,k], N0[j])
        p0[j,k] <- exp(lp[j,k])/(1+exp(lp[j,k]))
      lp[j,k] <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
      
      y1[j,k] ~ dbin(p1[j,k], N1[j])
        p1[j,k] <- exp(llp[j,k])/(1+exp(llp[j,k]))
      llp[j,k] <- b0[2] + b1[2]*trail[j] + b2[2]*cam1[j,k] + b3[2]*day1[j,k] + b4[2]*log[j]

    #probit(p0[j,k]) <- b0[1] + b1[1]*trail[j] + b2[1]*cam1[j,k] + b3[1]*day1[j,k] + b4[1]*log[j]
    } #k
  } #j

  #Derived parameters
  diff.slope <-  a1[1] - a1[2]
  diff.slope2 <-  a11[1] - a11[2]
  diff.elev <- a2[1] - a2[2]
  diff.elev2 <-  a3[1] - a3[2]
  diff.river <- a4[1] - a4[2]
  diff.log <-  a5[1] - a5[2]

} 
    ",fill = TRUE)
sink()

################################# BAYESIAN MODELING #############################
# Create the necessary arguments to run the bugs() command 
sp.data = list(y0=data0, y1=data1, nsite=nsite, nrep=dim(data0)[2], narea=narea, area=area,
               cam1=cam1, day1=day1, trail=trail, elev1=elev1, slope1=slope1, river1=river1, log=log) 

# Specify the parameters to be monitored
sp.params = list('a0', 'a1', 'a11', 'a2', 'a3', 'a4', 'a5', 'a6', 'b0', 'b1', 'b2', 'b3', 'b4', 'N0', 'N1', 'lambda0', 'lambda1',
	'diff.slope', 'diff.slope2', 'diff.elev', 'diff.elev2', 'diff.river', 'diff.log',
	'var.a7', 'var.a8', 'a7', 'a8')

# Specify the initial values
sp.inits = function() {
    list(a0=rnorm(2), a1=rnorm(2), a11=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(2), a6=rnorm(1),   
         b0=rnorm(2), b1=rnorm(2), b2=rnorm(2), b3=rnorm(2), b4=rnorm(2), a7=rnorm(narea), a8=rnorm(narea),
         sigma.a7=runif(1,0,3), sigma.a8=runif(1,0,3),
         N0 = as.vector(apply(data0,1,max, na.rm=T) + 5),   
         N1 = as.vector(apply(data1,1,max, na.rm=T) + 5)
    )           
}

# MCMC settings (choose one)
ni <- 300; nt <- 2; nb <- 100; nc <- 3 # quick test
ni <- 20000; nt <- 10; nb <- 10000; nc <- 3 # examining parameter values
ni <- 100000; nt <- 20; nb <- 50000; nc <- 3 # publication-quality run

# PORCUPINE
Sys.time()
out <- bugs(sp.data, sp.inits, sp.params, "covar3.txt", bugs.directory="c:/MyPrograms/WinBUGS14/",
	n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F)
save.image(file="workspace_porcupine.RData")
Sys.time()

#plot(out)
#out$summary  #Show mean, sd, CIs, second to last column is the Rhat, last column is the effect size

