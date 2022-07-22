


#==================================================================
# FULL MODEL
#==================================================================

# Libraries
library(ggplot2)
library(scales)
# Set working directory
setwd("~/Documents")


# Read in parameter values and ranges
params <- read.csv("params_enviro_120221.csv",header=F) 

# Read in experimental data points to fit to 
data <- read.csv("experimental_data_full_pens_WB_180520.csv",header=T) 

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
s_min <- params[21,3]
s_max <- params[22,3]
c_min <- params[23,3]
c_max <- params[24,3]
SE0_min <- 1
SE0_max <- 7
RE0_min <- 0
RE0_max <- 3
SC0_min <- 2
SC0_max <- 7
RC0_min <- 0
RC0_max <- 3

threshold_value <- function(th){
  # parameter values
  f <- f_min+th[1]*(f_max-f_min)
  Ke1 <- Ke_min+th[2]*(Ke_max-Ke_min)
  mue <- mue_min+th[3]*(mue_max-mue_min)
  beta_rs <- beta_rs_min+th[4]*(beta_rs_max-beta_rs_min)
  gse <- GSE_min+th[5]*(GSE_max-GSE_min)
  s <- s_min+th[6]*(s_max-s_min)
  c <- c_min+th[7]*(c_max-c_min)
  SE01 <- SE0_min+th[8]*(SE0_max-SE0_min)
  RE01 <- RE0_min+th[9]*(RE0_max-RE0_min)
  SC01 <- SC0_min+th[10]*(SC0_max-SC0_min)
  RC01 <- RC0_min+th[11]*(RC0_max-RC0_min)
  Ke <- (10^Ke1)
  SE <- (10^SE01)
  RE <- (10^RE01)
  SC <- (10^SC01)
  RC <- (10^RC01)
  
  # change scale from CFU/g to CFU/compartment
  E <- 10^7.542061
  A <- 30
  V <- 30
  Ke <- Ke*E
  RE <- RE*E
  SE <- SE*E
  SC <- SC*A*V*1000
  RC <- RC*A*V*1000
  c <- A*c/E

  # run model
  threshold <<- max(0,(RE + gse*(1-f)*RE*(1-(RE+SE)/Ke)+(beta_rs*(RE*SE)/(RE+SE))-mue*RE+s*RC-c*RE))/RE
  
}

param_names <- c("fitness cost", "carrying capacity", "decay rate", "horizontal transmission rate","growth rate","shedding rate","ingestion rate","density of non-resistant E. coli: enviro","density of resistant E. coli: enviro","density of non-resistant E. coli: cattle","density of resistant E. coli: cattle") 

# plot effect of individual parameters on R: plot 1
Percentiles_Uncertainty_Distribution <- numeric(0)
AMR_Reproductive_Rate <- numeric(0)
Variables <- numeric(0)
for(k in 1:11){
for(j in 1:1000){
theta <- runif(11,0,1.0)
thresholds <- numeric(11)
cutoffs <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
variables <- rep(param_names[k],11)
for(i in 1:11){
  theta[k] <- cutoffs[i]
  thresholds[i] <- threshold_value(theta)
}
Percentiles_Uncertainty_Distribution <- c(Percentiles_Uncertainty_Distribution,cutoffs)
AMR_Reproductive_Rate <- c(AMR_Reproductive_Rate,thresholds)
Variables <- c(Variables,variables)
}
}
data <- as.data.frame(cbind(Percentiles_Uncertainty_Distribution,AMR_Reproductive_Rate,Variables))
data <- data[order(Percentiles_Uncertainty_Distribution),]
data$Percentiles_Uncertainty_Distribution <- as.factor(data$Percentiles_Uncertainty_Distribution)
data$AMR_Reproductive_Rate <- as.numeric(data$AMR_Reproductive_Rate)
dim(data)
tail(data)
p <- ggplot(data, aes(x=Percentiles_Uncertainty_Distribution, y=AMR_Reproductive_Rate)) + 
  geom_violin()

# Plot effect of f on threshold, when all other parameters are medians
# create second dataset
Thresholds <- numeric(0)
Cutoffs <- numeric(0)
Variables <- numeric(0)
for(k in 1:11){
thresholds <- numeric(11)
theta <- rep(0.5,11)
for(i in 1:11){
theta[k] <- cutoffs[i]
thresholds[i] <- threshold_value(theta)
}
variables <- c(rep(param_names[k],11))

Thresholds <- c(Thresholds,thresholds)
Cutoffs <- c(Cutoffs,cutoffs)
Variables <- c(Variables,variables)
}
data2 <- as.data.frame(cbind(Cutoffs,Thresholds,Variables))
data2$Cutoffs <- as.factor(data2$Cutoffs)
data2$Thresholds <- as.numeric(data2$Thresholds)

p + geom_point(data=data2,aes(x=Cutoffs,y=Thresholds))+ facet_wrap(~ Variables, ncol = 4,nrow=3) + coord_cartesian(ylim = c(0,1.5),expand = TRUE,clip = "on"
) 

p + geom_point(data=data2,aes(x=Cutoffs,y=Thresholds))+ facet_wrap(~ Variables, ncol = 4,nrow=3)

