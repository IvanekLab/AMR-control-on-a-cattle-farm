


#==================================================================
# FIGURES
#==================================================================


# Colors
cols <- c("#000000","#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00")

# Packages needed
library(smfsb)

# Set working directory
setwd("~/Documents")

# Define functions for fitting ODEs
source(file = "Part3_equations_new_dates_130221.R") 
# Read in parameter values and ranges
params <- read.csv("params_enviro_120221.csv",header=F) 

# Read in experimental data points to fit to 
data <- read.csv("experimental_data_full_pens_WB_180520.csv",header=T) 

# read in temp data
temps <- read.csv("weather_data_cleaned_150921_WB.csv",header=T)

temps <- as.numeric(temps[(1:length(temps[,1])),1])
temps <- (temps-32)*5/9
length(temps[which(temps<(4))])

data_empty <- read.csv("experimental_data_enviro_WB_180520.csv",header=T) 
data_empty <- data_empty[c(1:5),]

# Define starting conditions
NE0 <- 10^(data$ML_generic[1]) 
RE0 <- 10^(data$ML_resistant[1])
SE0 <- NE0-RE0
NC0 <- 10^(data$ML_generic[6]) 
RC0 <- 10^(data$ML_resistant[6])
SC0 <- NC0-RC0

NE0_empty <- 10^data_empty$ML_generic[1]
RE0_empty <- 10^(data_empty$ML_resistant[1])
SE0_empty <- NE0_empty-RE0_empty

# Define experimental data points

t1_mean <- data$ML_generic[1]
t2_mean <- data$ML_generic[2]
t3_mean <- data$ML_generic[3]
t4_mean <- data$ML_generic[4]
t5_mean <- data$ML_generic[5]

t1_mean_cat <- data$ML_generic[6]
t2_mean_cat <- data$ML_generic[7]
t3_mean_cat <- data$ML_generic[8]
t4_mean_cat <- data$ML_generic[9]
t5_mean_cat <- data$ML_generic[10]

t1_Median_r <- data$ML_resistant[1]
t2_Median_r <- data$ML_resistant[2]
t3_Median_r <- data$ML_resistant[3]
t4_Median_r <- data$ML_resistant[4]
t5_Median_r <- data$ML_resistant[5]

t1_Median_r_cat <- data$ML_resistant[6]
t2_Median_r_cat <- data$ML_resistant[7]
t3_Median_r_cat <- data$ML_resistant[8]
t4_Median_r_cat <- data$ML_resistant[9]
t5_Median_r_cat <- data$ML_resistant[10]


t1_mean_empty <- data_empty$ML_generic[1]
t2_mean_empty <- data_empty$ML_generic[2]
t3_mean_empty <- data_empty$ML_generic[3]
t4_mean_empty <- data_empty$ML_generic[4]
t5_mean_empty <- data_empty$ML_generic[5]


t1_Median_r_empty <- data_empty$ML_resistant[1]
t2_Median_r_empty <- data_empty$ML_resistant[2]
t3_Median_r_empty <- data_empty$ML_resistant[3]
t4_Median_r_empty <- data_empty$ML_resistant[4]
t5_Median_r_empty <- data_empty$ML_resistant[5]


# define time-points
t1 <- data$day[1]
t2 <- data$day[2]
t3 <- data$day[3]
t4 <- data$day[4]
t5 <- data$day[5]

time_max <- (max(data$day))*24 # model runtime in days

param_names <- c("f", "Ke", "mue", "beta_rs","GSE","g","Kc","muc","s","c","I","n","Ke_empty","gse_empty","mue_empty","a") 

# Define fixed parameters
GSE_min <-params[1,3] # baseline growth rate in environment
GSE_max <- params[2,3]
f_min <- params[3,3] # fitness cost 
f_max <- params[4,3]
Ke_min <- params[5,3] # carrying capacity in the environment
Ke_max <- params[6,3]
mue_min <- params[7,3] # death rate in the environment
mue_max <- params[8,3]
beta_rs_min <- params[9,3] # horizantal transmission rate
beta_rs_max <- params[10,3]
g_min <- params[15,3]
g_max <- params[16,3]
Kc_min <- params[17,3]
Kc_max <- params[18,3]
muc_min <- params[19,3]
muc_max <- params[20,3]
s_min <- params[21,3]
s_max <- params[22,3]
c_min <- params[23,3]
c_max <- params[24,3]
I_min <- params[31,3]
I_max <- params[32,3]
n_min <- params[36,3] #remove carrying capacity change over time 
n_max <- params[37,3]

rmodel <- function(th){
  # parameter values
  f <- f_min+th[1]*(f_max-f_min)
  Ke <- Ke_min+th[2]*(Ke_max-Ke_min)
  mue <- mue_min+th[3]*(mue_max-mue_min)
  beta_rs <- beta_rs_min+th[4]*(beta_rs_max-beta_rs_min)
  gse <- GSE_min+th[5]*(GSE_max-GSE_min)
  g <- g_min+th[6]*(g_max-g_min)
  Kc <- Kc_min+th[7]*(Kc_max-Kc_min)
  muc <- muc_min+th[8]*(muc_max-muc_min)
  s <- s_min+th[9]*(s_max-s_min)
  c <- c_min+th[10]*(c_max-c_min)
  I <- I_min+th[11]*(I_max-I_min)
  n <- n_min+th[12]*(n_max-n_min)
  #I <- 0
  #n <- 0
  temperatures1 <- temps[1:(time_max+1)]
  gse_t <- temperatures1*I + gse
  #gse_t <- temperatures1*0 + gse # remove temperature from fitting procedure
  gse_t[which(temperatures1<4)] <- 0
  # run model
  resist_full.t <- seq(0,time_max,1)
  resist_full.init <- c(SE0,RE0,SC0,RC0) 
  pars <- c(f,beta_rs,g,Kc,muc,s,c,Ke,mue,n,gse_t)
  resist_full.sol <<- resist.dyn(resist_full.t, resist_full.init, pars)
  
  # label the output
  time <- resist_full.sol[,1]
  SE <- (resist_full.sol[,2]) # Take log10 of output
  RE <- (resist_full.sol[,3])
  NE <<- log10(SE+RE)
  RE <<- log10(RE)
  SC <- (resist_full.sol[,4]) # Take log10 of output
  RC <- (resist_full.sol[,5])
  NC <<- log10(SC+RC)
  RC <<- log10(RC)
  
  
  n <- 0
  s <- 0
  c <- 0
  Ke_empty <- 1.5+th[11]*(2.5-1.5) # change to having new parameters for empty pen: carrying capacity 1.5 to 2.5, baseline growth rate and mortality rate (same ranges)
  gse_empty <- GSE_min+th[12]*(GSE_max-GSE_min)# change vectors here when changing between basic and full model
  mue_empty <- mue_min+th[13]*(mue_max-mue_min)
  gse_t <- temperatures1*I + gse_empty
  gse_t[which(temperatures1<4)] <- 0
  
  #Ke_empty <- 2
  #mue_empty <- mue
  # run model
  resist_full.t <- seq(0,time_max,1)
  resist_full.init <- c(SE0_empty,RE0_empty,0,0) 
  pars <- c(f,beta_rs,g,Kc,muc,s,c,Ke_empty,mue_empty,n,gse_t)
  resist_full_empty.sol <<- resist.dyn(resist_full.t, resist_full.init, pars)
  SE_empty <- (resist_full_empty.sol[,2]) # Take log10 of output
  RE_empty <- (resist_full_empty.sol[,3])
  NE_empty <<- log10(SE_empty+RE_empty)
  RE_empty <<- log10(RE_empty)
}



#out <- read.csv("fit_params_full_model_WB_210122_basic.csv",header=T)# 20,000 iterations final model
out <- read.csv("fit_params_full_model_WB_210122.csv",header=T)# 20,000 iterations final model
out <- out[,c(2:17)]
out <- out[order(as.numeric(-out[,15])),]
head(out)
tail(out)


# print out fitted values (min and max)
f_min+out$f[1]*(f_max-f_min)
Ke_min+out$Ke[1]*(Ke_max-Ke_min)
mue_min+out$mue[1]*(mue_max-mue_min)
beta_rs_min+out$beta_rs[1]*(beta_rs_max-beta_rs_min)
gse <- GSE_min+out$GSE[1]*(GSE_max-GSE_min)
g_min+out$g[1]*(g_max-g_min)
Kc_min+out$Kc[1]*(Kc_max-Kc_min)
muc_min+out$muc[1]*(muc_max-muc_min)
s_min+out$s[1]*(s_max-s_min)
c_min+out$c[1]*(c_max-c_min)
I <- I_min+out$I[1]*(I_max-I_min)
n_min+out$n[1]*(n_max-n_min)
gse+I*4

f_min+min(out$f)*(f_max-f_min)
Ke_min+min(out$Ke)*(Ke_max-Ke_min)
mue_min+min(out$mue)*(mue_max-mue_min)
beta_rs_min+min(out$beta_rs)*(beta_rs_max-beta_rs_min)
gse <- GSE_min+min(out$GSE)*(GSE_max-GSE_min)
g_min+min(out$g)*(g_max-g_min)
Kc_min+min(out$Kc)*(Kc_max-Kc_min)
muc_min+min(out$muc)*(muc_max-muc_min)
s_min+min(out$s)*(s_max-s_min)
c_min+min(out$c)*(c_max-c_min)
I <- I_min+min(out$I)*(I_max-I_min)
n_min+min(out$n)*(n_max-n_min)
gse+I*4

f_min+max(out$f)*(f_max-f_min)
Ke_min+max(out$Ke)*(Ke_max-Ke_min)
mue_min+max(out$mue)*(mue_max-mue_min)
beta_rs_min+max(out$beta_rs)*(beta_rs_max-beta_rs_min)
gse <- GSE_min+max(out$GSE)*(GSE_max-GSE_min)
g_min+max(out$g)*(g_max-g_min)
Kc_min+max(out$Kc)*(Kc_max-Kc_min)
muc_min+max(out$muc)*(muc_max-muc_min)
s_min+max(out$s)*(s_max-s_min)
c_min+max(out$c)*(c_max-c_min)
I <- I_min+max(out$I)*(I_max-I_min)
n_min+max(out$n)*(n_max-n_min)
gse+I*4

interventions <- out[1,c(1:12)]
interventions[2] <- 0

interventions1 <- out[1,c(1:12)]
interventions1[2] <- 0

# Figure of model fit 
png("Fig2.png",400,400)  

par(mfrow=c(1,1))

plot(0,xlim=c(0,(126*24)),ylim=c(-4,10),type='n',xlab="Hours since cattle introduced to pen",ylab="Log10 CFU/g",main=substitute(paste("Tetracycline-resistant ", italic("E. coli"))))

rmodel(as.numeric(out[1,c(1:16)])) # change vector here
points(c(1:length(RE)),RE,col=cols[1],pch=19,cex=0.25)
points(c(1:length(RC)),RC,col=cols[2],pch=19,cex=0.25)

points(data$day[1:5]*24,data$ML_resistant[1:5],col=cols[1])
points(data$day[6:10]*24,data$ML_resistant[6:10],col=cols[2])

arrows(data$day[1:5]*24,data$Lower_CI_res[1:5],data$day[1:5]*24,data$Upper_CI_res[1:5],col=cols[1],angle=90,code=3)
arrows(data$day[6:10]*24,data$Lower_CI_res[6:10],data$day[6:10]*24,data$Upper_CI_res[6:10],col=cols[2],angle=90,code=3)
legend(1000, 0, legend=c("model predictions:cattle", "model predictions:environment", "trial data:cattle","trial data:environment"),
       col=c(cols[2],cols[1]), lty =c(1,1,NA,NA),pch=c(NA,NA,1,1), cex=1)
dev.off() 
plot(0,xlim=c(0,(126*24)),ylim=c(-4,10),type='n',xlab="hours",ylab="Log10 CFU/g",main="Total Generic E. coli")

points(c(1:length(NE)),NE,col=cols[1],pch=19,cex=0.25)
points(c(1:length(NC)),NC,col=cols[2],pch=19,cex=0.25)

points(data$day[1:5]*24,data$ML_generic[1:5],col=cols[1])
points(data$day[6:10]*24,data$ML_generic[6:10],col=cols[2])

arrows(data$day[1:5]*24,data$Lower_CI_generic[1:5],data$day[1:5]*24,data$Upper_CI_generic[1:5],col=cols[1],angle=90,code=3)
arrows(data$day[6:10]*24,data$Lower_CI_generic[6:10],data$day[6:10]*24,data$Upper_CI_generic[6:10],col=cols[2],angle=90,code=3)

# Look at empty pen predictions
# Read in experimental data points to fit to 
png("Fig3.png",400,400) 
plot(0,xlim=c(0,(126*24)),ylim=c(-4,10),type='n',xlab="Hours since cattle introduced to pen",ylab="Log10 CFU/g",main=substitute(paste("Tetracycline-resistant ", italic("E. coli"))))

points(c(1:length(RE_empty)),RE_empty,col=cols[1],pch=19,cex=0.25)

points(data_empty$day[1:5]*24,data_empty$ML_resistant[1:5],col=cols[1]) # change to enviro data
arrows(data_empty$day[1:5]*24,data_empty$Lower_CI_res[1:5],data_empty$day[1:5]*24,data_empty$Upper_CI_res[1:5],col=cols[1],angle=90,code=3)
legend(1100, 0.75, legend=c( "Model prediction", "Trial data from empty pens"),
       col=c(cols[1],cols[1]), lty =c(1,NA), pch=c(NA,1), cex=0.85)

dev.off() 
plot(0,xlim=c(0,(126*24)),ylim=c(-4,10),type='n',xlab="hours",ylab="Log10 CFU/g",main="Total Generic E. coli")

points(c(1:length(NE_empty)),NE_empty,col=cols[1],pch=19,cex=0.25)

points(data_empty$day[1:5]*24,data_empty$ML_generic[1:5],col=cols[1])
arrows(data_empty$day[1:5]*24,data_empty$Lower_CI_generic[1:5],data_empty$day[1:5]*24,data_empty$Upper_CI_generic[1:5],col=cols[1],angle=90,code=3)

legend(1100, 0.75, legend=c( "Model prediction", "Trial data from empty pens"),
       col=c(cols[1],cols[1]), lty =c(1,NA), pch=c(NA,1), cex=0.85)
# Figure of interventions
par(mfrow=c(1,1))
png("Fig4.png",400,400,pointsize = 12)  ## open a pdf file for plotting

plot(0,xlim=c(0,(126*24)),ylim=c(-4,12),type='n',xlab="Hours since cattle introduced to pen",ylab="Log10 CFU/g",main=substitute(paste("Tetracycline-resistant ", italic("E. coli"))))

rmodel(as.numeric(out[1,c(1:13)]))
lines(c(1:length(RE)),RE,col=cols[1],pch=19,cex=0.25,lty=1)
lines(c(1:length(RC)),RC,col=cols[2],pch=19,cex=0.25,lty=1)


rmodel(as.numeric(interventions))
lines(c(1:length(RE)),RE,col=cols[1],pch=19,cex=0.25,lty=2)
lines(c(1:length(RC)),RC,col=cols[2],pch=19,cex=0.25,lty=2)

Ke_min <- 1 # carrying capacity in the environment

rmodel(as.numeric(interventions1))
lines(c(1:length(RE)),RE,col=cols[1],pch=19,cex=0.25,lty=3)
lines(c(1:length(RC)),RC,col=cols[2],pch=19,cex=0.25,lty=3)
Ke_min <- params[5,3] # carrying capacity in the environment

legend(5, 12, legend=c( "Cattle: baseline (no scraping)", "Environment: baseline","Cattle: increased pen scraping","Environment: increased pen scraping","Cattle: maximum pen scraping","Environment: maximum pen scraping"),
       col=c(cols[2],cols[1]), lty =c(1,1,2,2,3,3), cex=1)
dev.off()
plot(0,xlim=c(0,(126*24)),ylim=c(-2,7.5),type='n',xlab="hours",ylab="Log10 CFU/g",main="Total Generic E. coli per gram")

rmodel(as.numeric(out[1,c(1:13)]))
lines(c(1:length(NE)),NE,col=cols[1],pch=19,cex=0.25,lty=1)
lines(c(1:length(NC)),NC,col=cols[2],pch=19,cex=0.25,lty=1)

rmodel(as.numeric(interventions))
lines(c(1:length(NE)),NE,col=cols[1],pch=19,cex=0.25,lty=2)
lines(c(1:length(NC)),NC,col=cols[2],pch=19,cex=0.25,lty=2)
Ke_min <- 1 # carrying capacity in the environment

rmodel(as.numeric(interventions1))
lines(c(1:length(NE)),NE,col=cols[1],pch=19,cex=0.25,lty=3)
lines(c(1:length(NC)),NC,col=cols[2],pch=19,cex=0.25,lty=3)
Ke_min <- params[5,3] # carrying capacity in the environment


dev.off() 
# Figure of relationship between temperature and g
png("Fig5.png",400,400,pointsize = 12)  
par(mfrow=c(1,1))

tempC <- c(4:35)
gs <- GSE_min + out[1,5]*(GSE_max-GSE_min)
i <- I_min+ out[1,11]*(I_max-I_min)
plot(0,xlim=c(0,35),ylim=c(-0.01,3),type='n',xlab="Ambient temperature in Celsius",ylab="Predicted Specific Maximum Growth Rate of E. coli ")
lines(tempC, i*tempC+gs)
# Add points from Gautam: E. coli O157 in cattle manure
gs <- 0.0032
i <- -0.0163
lines(tempC, gs*tempC+i, col=cols[2])
legend(0, 2.5, legend=c("Best fit values- generic E. coli","Observed values Gautam et al. 2011 - E. coli O157:H7"),
       col=c(cols[1],cols[2]), lty =c(1,1), cex=0.85)
dev.off() 
# Figure of predictions of fitted model at different constant temperatures 
par(mfrow=c(1,1))

png("Fig6.png",400,400)  ## open a pdf file for plotting

plot(0,xlim=c(0,(126*24)),ylim=c(-10,10),type='n',xlab="hours",ylab="Log10 CFU/g",main="Tetracycline-resistant E. coli")
temps <- rep(30,3384)
rmodel(as.numeric(out[1,c(1:13)]))
lines(c(1:length(RE)),RE,col=cols[1],lty=1,cex=0.25)
lines(c(1:length(RC)),RC,col=cols[2],lty=1,cex=0.25)

temps <- rep(0,3384)
rmodel(as.numeric(out[1,c(1:13)]))
lines(c(1:length(RE)),RE,col=cols[1],lty=2,cex=0.25)
lines(c(1:length(RC)),RC,col=cols[2],lty=2,cex=0.25)

legend(300, 4, legend=c("Cattle ingesta: ambient temperature 30 degrees C", "Pen environment: 30 degrees C","Cattle ingesta: ambient temperature 0 degrees C","Pen environment: 0 degrees C"),
       col=c(cols[2],cols[1]), lty =c(1,1,2,2), cex=0.85)

dev.off() 

plot(0,xlim=c(0,(126*24)),ylim=c(-10,10),type='n',xlab="hours",ylab="Log10 CFU/g",main="Total Generic E. coli")
temps <- rep(30,3384)
rmodel(as.numeric(out[1,c(1:15)]))
lines(c(1:length(NE)),NE,col=cols[1],lty=1,cex=0.25)
lines(c(1:length(NC)),NC,col=cols[2],lty=1,cex=0.25)

temps <- rep(0,3384)
rmodel(as.numeric(out[1,c(1:12)]))
lines(c(1:length(NE)),NE,col=cols[1],lty=2,cex=0.25)
lines(c(1:length(NC)),NC,col=cols[2],lty=2,cex=0.25)

# Plot temperatures
par(mfrow=c(1,1))
temps <- read.csv("weather_data_cleaned_150921_WB.csv",header=T)

temps <- as.numeric(temps[(1:length(temps[,1])),1])
temps <- (temps-32)*5/9
plot(1:3024,temps[1:3024])

