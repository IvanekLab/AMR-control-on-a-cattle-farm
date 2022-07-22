


#==================================================================
# FULL MODEL
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
length(temps)
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
  
  
  n <- th[12]*0
  s <- 0
  c <- 0
  Ke_empty <- 1.5+th[13]*(2.5-1.5) # change to having new parameters for empty pen: carrying capacity 1.5 to 2.5, baseline growth rate and mortality rate (same ranges)
  gse_empty <- GSE_min+th[14]*(GSE_max-GSE_min)
  mue_empty <- mue_min+th[15]*(mue_max-mue_min)
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

# Difference function - use model function

rdist <- function(th){
  rmodel(th)
  sum(sqrt((RE[(t2*24+1)]-t2_Median_r)^2)+
      sqrt((RE[(t3*24+1)]-t3_Median_r)^2)+
      sqrt((RE[(t4*24+1)]-t4_Median_r)^2)+
      sqrt((RE[(t5*24+1)]-t5_Median_r)^2)+
        sqrt((NE[(t2*24+1)]-t2_mean)^2)+
      sqrt((NE[(t3*24+1)]-t3_mean)^2)+
      sqrt((NE[(t4*24+1)]-t4_mean)^2)+
      sqrt((NE[(t5*24+1)]-t5_mean)^2)+
        sqrt((RC[(t2*24+1)]-t2_Median_r_cat)^2)+
      sqrt((RC[(t3*24+1)]-t3_Median_r_cat)^2)+
      sqrt((RC[(t4*24+1)]-t4_Median_r_cat)^2)+
      sqrt((RC[(t5*24+1)]-t5_Median_r_cat)^2)+
        sqrt((NC[(t2*24+1)]-t2_mean_cat)^2)+
      sqrt((NC[(t3*24+1)]-t3_mean_cat)^2)+
      sqrt((NC[(t4*24+1)]-t4_mean_cat)^2)+
      sqrt((NC[(t5*24+1)]-t5_mean_cat)^2)+
        sqrt((RE_empty[(t2*24+1)]-t2_Median_r_empty)^2)+
    sqrt((RE_empty[(t3*24+1)]-t3_Median_r_empty)^2)+
    sqrt((RE_empty[(t4*24+1)]-t4_Median_r_empty)^2)+
    sqrt((RE_empty[(t5*24+1)]-t5_Median_r_empty)^2)+
    sqrt((NE_empty[(t2*24+1)]-t2_mean_empty)^2)+
    sqrt((NE_empty[(t3*24+1)]-t3_mean_empty)^2)+
    sqrt((NE_empty[(t4*24+1)]-t4_mean_empty)^2)+
    sqrt((NE_empty[(t5*24+1)]-t5_mean_empty)^2),na.rm=T)
  
  
}

rprior <- function() { 
  
  # parameter values
  runif(15,0,1)
}


dprior <- function(x, ...) { dunif(x[1], 0,1, ...) + 
    dunif(x[2], 0,1, ...) + dunif(x[3], 0,1, ...) + dunif(x[4], 0,1, ...)+ dunif(x[5], 0,1, ...)+ dunif(x[6], 0,1, ...)+ dunif(x[7], 0,1, ...)+ dunif(x[8], 0,1, ...)+ dunif(x[9], 0,1, ...)+ dunif(x[10], 0,1, ...)+ dunif(x[11], 0,1, ...)+ dunif(x[12], 0,1, ...)+ dunif(x[13], 0,1, ...)+ dunif(x[14], 0,1, ...)+ dunif(x[15], 0,1, ...)}

rperturb <- function(th){th + rnorm(15, 0, 0.005)}


dperturb <- function(thNew, thOld, ...){sum(dnorm(thNew, thOld, 0.005, ...))} 

# Model settings: 
total_iterations <- 100 
start <-proc.time()
out = abcSmc(total_iterations, rprior, dprior, rdist, rperturb,
             dperturb, verb=T, steps=10, factor=10)

a <- numeric(total_iterations)
for(i in 1:total_iterations){
  a[i] <- rdist(out[i,])
}
out <- cbind(out,a)
colnames(out) <- c(param_names)


# save output
#write.csv(out,paste("output.csv"))
#print(proc.time()-start)
