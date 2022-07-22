
#==================================================================
# FULL MODEL: ENVIRONMENT & ANIMALS
#==================================================================


#==================================================================
#Model Equations 
#==================================================================
resist.dyn <- function(time,var,pars) {
  # Input parameters 
  #gse <- pars[1] # baseline growth rate in environment
  f <- pars[1] # fitness cost 
  beta_rs <- pars[2] # horizantal transfer rate in the environment & in cattle
  gsc <- pars[3] # baseline growth rate in cattle
  Kc <- pars[4] # carrying capacity in cattle (total cattle compartment)
  muc <- pars[5] # death rate in the cattle
  s <- pars[6] # rate at which bacteria are shed from cattle to environment (via feces)
  c <- pars[7] # rate at which bacteria are spread from environment to cattle (via consumption)
  Ke <-pars[8]
  mue <-pars[9]
  n <-pars[10]
  gse <- pars[11:3036]

  # Calculated parameters
  grc <- gsc*(1-f)
  
  # States
  SE <- var[1] # Non-resistant bacteria in the environment
  RE <- var[2] # Resistant bacteria in the environment
  SC <- var[3] # Non-resistant bacteria in the cattle
  RC <- var[4] # Resistant bacteria in the cattle
  
  #Create vectors
  SE_final <- numeric(max(time)+1)
  RE_final <- numeric(max(time)+1)
  SC_final <- numeric(max(time)+1)
  RC_final <- numeric(max(time)+1)
  
  SE_final[1] <- (SE)
  RE_final[1] <- (RE)
  SC_final[1] <- (SC)
  RC_final[1] <- (RC)
  
  Ke <- (10^Ke)
  Kc <- (10^Kc)
  
  
  E <- 10^7.542061
  A <- 30
  V <- 30
  ke <- Ke*E
  Kc <- Kc*A*V*1000
  RE <- RE*E
  SE <- SE*E
  SC <- SC*A*V*1000
  RC <- RC*A*V*1000
  c <- A*c/E

  
  # Calculate the next time-step
  for(i in 1:max(time)){
  SE1 <- max(0,SE + gse[i]*SE*(1-(RE+SE)/(ke+n*i))-mue*SE) # growth and death
  RE1 <- max(0,RE + gse[i]*(1-f)*RE*(1-(RE+SE)/(ke+n*i))-mue*RE)
  SC1 <- max(0,SC + gsc*SC*(1-(RC+SC)/Kc)-muc*SC)
  RC1 <- max(0,RC + grc*RC*(1-(RC+SC)/Kc)-muc*RC)
  
  SE2 <- SE1-max((beta_rs*(RE1*SE1)/(RE1+SE1)),0,na.rm=T) 
  RE2 <- RE1+max((beta_rs*(RE1*SE1)/(RE1+SE1)),0,na.rm=T)# horizantal transmission
  SC2 <- SC1-max((beta_rs*(RC1*SC1)/(RC1+SC1)),0,na.rm=T)
  RC2 <- RC1+max((beta_rs*(RC1*SC1)/(RC1+SC1)),0,na.rm=T)
  
  SE3 <- SE2+s*SC2 # shedding
  RE3 <- RE2+s*RC2
  SC3 <- SC2-s*SC2
  RC3 <- RC2-s*RC2
  
  SE <- SE3-c*SE3# ingestion
  RE <- RE3-c*RE3
  SC <- SC3+c*SE3
  RC <- RC3+c*RE3
  
  SE_final[i+1] <- (SE/E)
  RE_final[i+1] <- (RE/E)
  SC_final[i+1] <- (SC/(A*V*1000))
  RC_final[i+1] <- (RC/(A*V*1000))
  }
  # Last instruction: return a list 
  return((cbind(time,SE_final,RE_final,SC_final,RC_final)))
}



