### microPRIME functions ###

# Last updated: 12th March 2021

# 'assign' function removed from trend functions on 06/01/21

# The purpose of this R script is to store all of the functions that are used in the microPRIME
# model. Each of the parameters used in the models are provided as point estimates here. There
# is the opportunity to allow any of the parameters to be varying or unvarying between iterations.
# This is set in the microsim_multiple_popstill_parallel.R script, just after the loop starts.
# The varying parameters go forward to the emulation phase. The unvarying parameters are 'hard-wired'
# into the model process

# 1. Parameter point estimates
# 2a. Relative risk functions
# 2b. Risk factor trend functions
# 2c. State transition functions

##############################################################################
#                                                                            #
# 1. Parameter point estimates                                               #
#                                                                            #
##############################################################################

# All variables below have ".pe" signifying a point estimate. When not selected as a varying
# parameter, these values will be selected to be used in the microPRIME functions.

minrisk.bmi.pe <- 15        # This sets the minimum risk level for BMI, where by definition RR is 1.
sbplow.pe      <- 90        # This sets the minimum risk level for blood pressure, where by definition RR is 1.
tcholow.pe     <- 2         # This sets the minimum risk level for total cholesterol, where by definition RR is 1.
combrr.max.pe  <- 2500      # This sets the maximum combined relative risk for any individual. Effectively, this 

# See accompanying Excel files for determining shape of transition functions, relative risks and risk factor
# trends. 
alpha.30dcf.m.pe          <- -2.8473 # Case-fatality
beta.30dcf.m.pe           <-  0.0285 
alpha.30dcf.f.pe          <- -2.7306
beta.30dcf.f.pe           <-  0.0267 
alpha.mort.m.pe           <- 32.1538 # See demo.R (coefs2)
beta1.mort.m.pe           <-  0.0817 # Exponential background mortality curve
beta2.mort.m.pe           <- -0.0208
alpha.mort.f.pe           <- 24.7410
beta1.mort.f.pe           <-  0.0837
beta2.mort.f.pe           <- -0.0174
treatment.inc.mi.m.pe     <-  0.0000 #0.0200 # 2% annual decline in incidence
treatment.cf.mi.m.pe      <-  0.0500 # 5% annual decline in case fatality
treatment.inc.mi.f.pe     <-  0.0000 #0.0200 # NB: Set to zero for runs without treatment effect on incidence
treatment.cf.mi.f.pe      <-  0.0500
bmi.m.r1525.pe            <-  1.2700 # BMI RR
bmi.m.r2550.pe            <-  1.4200
bmi.f.r1525.pe            <-  1.0100
bmi.f.r2550.pe            <-  1.3500
sbp.alpha.m.pe            <-  0.2470 # SBP RR
sbp.beta.m.pe             <-  0.0050 
sbp.alpha.f.pe            <-  0.1650
sbp.beta.f.pe             <-  0.0054
tcho.alpha.m.pe           <- -1.8938 # Cholesterol RR
tcho.beta.m.pe            <-  0.0407 
tcho.alpha.f.pe           <- -3.5567
tcho.beta.f.pe            <-  0.0706
smok.alpha.m.pe           <-  4.1338 # Smoking RR
smok.beta.m.pe            <- -0.0310 
smok.alpha.f.pe           <-  5.2915
smok.beta.f.pe            <- -0.0423
diab.m.pe                 <-  1.8500 # Diabetes RR
diab.f.pe                 <-  2.6300
prev.alpha.m.pe           <-  6.0431 # Previous MI RR 
prev.alpha.f.pe           <-  7.3793
prev.beta.age.m.pe        <- -0.0402
prev.beta.age.f.pe        <- -0.0656
prev.beta.time.m.pe       <- -0.4121
prev.beta.time.f.pe       <- -0.4404
sbp.trend.age.m.pe        <-  1.3805 # SBP trend
sbp.trend.age.f.pe        <-  1.8740
sbp.trend.age.quad.m.pe   <- -0.0080
sbp.trend.age.quad.f.pe   <- -0.0092
sbp.trend.year.m.pe       <- -0.9890
sbp.trend.year.f.pe       <- -1.1330
sbp.trend.year.quad.m.pe  <-  0.0116
sbp.trend.year.quad.f.pe  <-  0.0089
tcho.trend.age.m.pe       <-  0.1212 # Cholesterol trend
tcho.trend.age.f.pe       <-  0.0990
tcho.trend.age.quad.m.pe  <- -0.0011
tcho.trend.age.quad.f.pe  <- -0.0007
tcho.trend.year.m.pe      <- -0.0339
tcho.trend.year.f.pe      <- -0.0473
tcho.trend.year.quad.f.pe <-  0.0006
smok.trend.alpha.m.pe     <- -0.3142 # Smoking trend 
smok.trend.age.m.pe       <-  0.0121
smok.trend.age.quad.m.pe  <- -0.0004
smok.trend.year.m.pe      <- -0.0285
smok.trend.alpha.f.pe     <- -0.5873
smok.trend.age.f.pe       <-  0.0155
smok.trend.age.quad.f.pe  <- -0.0004
smok.trend.year.f.pe      <- -0.0205
smok.trend.year.quad.f.pe <- -0.0005
diab.trend.alpha.m.pe     <- -9.9950 # Diabetes trend
diab.trend.age.m.pe       <-  0.1804
diab.trend.age.quad.m.pe  <- -0.0011
diab.trend.alpha.f.pe     <- -8.8070
diab.trend.age.f.pe       <-  0.1323
diab.trend.age.quad.f.pe  <- -0.0007
diab.trend.year.m.pe      <-  0.0147
diab.trend.year.quad.m.pe <-  0.0012
diab.trend.year.f.pe      <-  0.0483

## Best case / worst case unvarying parameters for BMI trends
bmi.meantrend.const.pe     <-  2.972980 # Best case, based on linear polynomial models
bmi.meantrend.sex.pe       <-  0.014412
bmi.meantrend.age.pe       <-  0.011427
bmi.meantrend.year.pe      <-  0.002521
bmi.meantrend.age.quad.pe  <- -0.000095
bmi.meantrend.year.quad.pe <- -0.000121
bmi.sdtrend.const.pe       <- -3.806541
bmi.sdtrend.sex.pe         <- -0.023976
bmi.sdtrend.mean.pe        <-  2.452630
bmi.sdtrend.mean.quad.pe   <- -0.373192
bmi.sdtrend.age.pe         <- -0.000382
bmi.sdtrend.year.pe        <-  0.002788
bmi.sdtrend.age.quad.pe    <- -0.000003
bmi.sdtrend.year.quad.pe   <- -0.000041
bmi.sdtrend.sexage.pe      <- -0.000937
bmi.sdtrend.sexagesq.pe    <-  0.000010
bmi.sdtrend.ageyear.pe     <- -0.000035
bmi.sdtrend.ageyearsq.pe   <-  0.000001
bmi.meannl.a.pe            <-  3.005400 # Worst case, based on non-linear models
bmi.meannl.b.pe            <-  0.024797
bmi.meannl.year.pe         <-  0.090790
bmi.meannl.sex.pe          <- -0.002427
bmi.meannl.age.pe          <-  0.011050
bmi.meannl.age.quad.pe     <- -0.000092
bmi.meannl.sexage.pe       <-  0.000787
bmi.meannl.sexagesq.pe     <- -0.000007
bmi.sdnl.a.pe              <-  0.176627
bmi.sdnl.b.pe              <- -0.018413
bmi.sdnl.year.pe           <- -0.048294
bmi.sdnl.sex.pe            <- -0.023324
bmi.sdnl.age.pe            <-  0.000052
bmi.sdnl.age.quad.pe       <- -0.000006
bmi.sdnl.sexage.pe         <- -0.000973
bmi.sdnl.sexagesq.pe       <-  0.000010
bmi.sdnl.mean.pe           <-  0.004284



##############################################################################
#                                                                            #
# 2a. Relative risk functions                                                #
#                                                                            #
##############################################################################

# Function that selects the varying or point estimate parameters
take <- function(z) {
  df$modelparameter[df$name==z]
}


# BMI
rrbmi.m <- function (bmi) { # Function for men
  if (bmi<take("minrisk.bmi")) {rr <- 1}     
  else if (bmi<25) {rr <- take("bmi.m.r1525")^((bmi-take("minrisk.bmi"))/5)}
  else {rr <- take("bmi.m.r1525")^((25-take("minrisk.bmi"))/5)*take("bmi.m.r2550")^((bmi-25)/5)} 
  return(rr)
}

rrbmi.f <- function (bmi) { # Function for women
  if (bmi<take("minrisk.bmi")) {rr <- 1}     
  else if (bmi<25) {rr <- take("bmi.f.r1525")^((bmi-take("minrisk.bmi"))/5)}
  else {rr <- take("bmi.f.r1525")^((25-take("minrisk.bmi"))/5)*take("bmi.f.r2550")^((bmi-25)/5)} 
  return(rr)
}

rrbmi <- function (sex, bmi) {
  if (sex==2) {rr <- rrbmi.m(bmi)}      # Men
  else if (sex==1) {rr <- rrbmi.f(bmi)} # Women
  else {rr <- c("Error")}
  return(rr)
}

# SBP
rrsbp <- function (age, sex, sbp) {
  if (sbp < take("sbplow")) {rr <- 1}
  else if (sex==1) {rr <- max(0,(take("sbp.alpha.f")+take("sbp.beta.f")*age)^((take("sbplow")-sbp)/20))}
  else if (sex==2) {rr <- max(0,(take("sbp.alpha.m")+take("sbp.beta.m")*age)^((take("sbplow")-sbp)/20))}
  else {rr <- c("Error")}
  return(rr)
}

# Total cholesterol
sigmoid.fun <- function(a,b,age) {
  1/(1+exp(-(a+b*age)))
}

rrtcho <- function (age, sex, tcho) {
  if (tcho<take("tcholow")) {rr <- 1}
  else if (sex==1) {rr <- max(0,sigmoid.fun(take("tcho.alpha.f"),take("tcho.beta.f"),age)^((take("tcholow")-tcho)/1))}
  else if (sex==2) {rr <- max(0,sigmoid.fun(take("tcho.alpha.m"),take("tcho.beta.m"),age)^((take("tcholow")-tcho)/1))}
  else {rr <- c("Error")}
  return(rr)
}

# Smoking
rrsmok <- function (age, sex, smok) {
  if (smok==0) {rr <- 1}
  else if (sex==1) {rr <- max(0,take("smok.alpha.f")+take("smok.beta.f")*age)}
  else if (sex==2) {rr <- max(0,take("smok.alpha.m")+take("smok.beta.m")*age)}
  else {rr <- c("Error")}
  return(rr)
}

# Diabetes
rrdiab <- function(sex, diab) {
  if(diab==0) {rr <- 1}
  else if (sex==1) {rr <- take("diab.f")}
  else if (sex==2) {rr <- take("diab.m")}
  else {rr <- c("Error")}
  return(rr)
}

# Previous MI relative risk
rrprev <- function(age, sex, time.mi) {
  if (time.mi==0) {rr <- 1}
  else if (sex==1) {rr <- max(1,exp(take("prev.alpha.f")+take("prev.beta.age.f")*age+take("prev.beta.time.f")*time.mi))}
  else if (sex==2) {rr <- max(1,exp(take("prev.alpha.m")+take("prev.beta.age.m")*age+take("prev.beta.time.m")*time.mi))}
  else {rr <- c("Error")}
  return(rr)
}

# Function for combining RRs (without previous MIs as not available in baseline rate cohort)
rr.fun <- function(sex,age,bmi,sbp,tcho,smok,diab) {
  min(take("combrr.max"),
      rrbmi(sex,bmi)*rrsbp(age,sex,sbp)*rrtcho(age,sex,tcho)*rrsmok(age,sex,smok)*rrdiab(sex,diab))
}

##############################################################################
#                                                                            #
# 2b. Risk factor trend functions                                            #
#                                                                            #
##############################################################################

## Note that 1991 is hard-wired into the trend equations for sbp, tcho, smok and diab,
## and 2003 is hard-wired into trend equations for BMI.
## BMI trend equations described fully in Cobiac L, Scarborough P. Modelling future trajectories
## of obesity and body mass index in England (in progress).
## Continuous risk factors for sbp and tcho work by adding the yearly change to current 
## risk factor for each individual.
## For BMI, the z score for each individual (which remains constant) is translated into an 
## actual BMI value

# BMI
z.fun <- function(bmi, meanlnbmi, sdlnbmi) { # This transforms the initial bmi into a z score for projection
  (log(bmi) - meanlnbmi)/sdlnbmi
}

backz <- function(z, mean, sd) { # This transforms the z score to an actual value
  exp(mean + sd*z)
}

if (bestworst == "best") {
  mean.bmi.trend <- function(age, sex, year) {
    if (sex == 1) { # Note that for Cobiac et al parameters, male = 1, female = 0
      sum(take("bmi.meantrend.const"),
          take("bmi.meantrend.sex")*sex,
          take("bmi.trend.age.m")*age,
          take("bmi.meantrend.age.quad")*age*age,
          take("bmi.trend.year.m")*year,
          take("bmi.meantrend.year.quad")*year*year)
    } else if (sex==0) {
      sum(take("bmi.meantrend.const"),
          take("bmi.meantrend.sex")*sex,
          take("bmi.trend.age.f")*age,
          take("bmi.meantrend.age.quad")*age*age,
          take("bmi.trend.year.f")*year,
          take("bmi.meantrend.year.quad")*year*year)
    } else { print("Error in sex in bmi trend function") }
  }
  sd.bmi.trend <- function(age, sex, meanbmi, year) {
    if (sex==1) { # Note that for Cobiac et al parameters, male = 1, female = 0
      sum(take("bmi.sdtrend.const"),
          take("bmi.sdtrend.sex")*sex,
          take("bmi.sdtrend.age.m")*age,
          take("bmi.sdtrend.year.m")*year,
          take("bmi.sdtrend.age.quad")*age*age,
          take("bmi.sdtrend.year.quad")*year*year,
          take("bmi.sdtrend.sexage")*sex*age,
          take("bmi.sdtrend.sexagesq")*sex*age*age,
          take("bmi.sdtrend.ageyear")*age*year,
          take("bmi.sdtrend.ageyearsq")*age*year*year,
          take("bmi.sdtrend.mean")*meanbmi,
          take("bmi.sdtrend.mean.quad")*meanbmi*meanbmi)
    } else if (sex == 0) {
      sum(take("bmi.sdtrend.const"),
          take("bmi.sdtrend.sex")*sex,
          take("bmi.sdtrend.age.f")*age,
          take("bmi.sdtrend.year.f")*year,
          take("bmi.sdtrend.age.quad")*age*age,
          take("bmi.sdtrend.year.quad")*year*year,
          take("bmi.sdtrend.sexage")*sex*age,
          take("bmi.sdtrend.sexagesq")*sex*age*age,
          take("bmi.sdtrend.ageyear")*age*year,
          take("bmi.sdtrend.ageyearsq")*age*year*year,
          take("bmi.sdtrend.mean")*meanbmi,
          take("bmi.sdtrend.mean.quad")*meanbmi*meanbmi)
    } else { print("Error in sex in bmi trend function") }
  }
} else if (bestworst == "worst") {
  mean.bmi.trend <- function(age, year, sex) {
    if (sex == 1) { # Note that for Cobiac et al parameters, male = 1, female = 0
      sum(take("bmi.meannl.a"),
          -take("bmi.meannl.b")*exp(-take("bmi.trend.year.m")*year),
          take("bmi.meannl.sex")*sex,
          take("bmi.trend.age.m")*age,
          take("bmi.meannl.age.quad")*age*age,
          take("bmi.meannl.sexage")*sex*age,
          take("bmi.meannl.sexagesq")*sex*age*age)
    } else if (sex == 0) {
      sum(take("bmi.meannl.a"),
          -take("bmi.meannl.b")*exp(-take("bmi.trend.year.f")*year),
          take("bmi.meannl.sex")*sex,
          take("bmi.trend.age.f")*age,
          take("bmi.meannl.age.quad")*age*age,
          take("bmi.meannl.sexage")*sex*age,
          take("bmi.meannl.sexagesq")*sex*age*age)
    } else { print("Error in sex in bmi trend function") }
  }
  sd.bmi.trend <- function(age, year, sex, meanbmi) {
    if (sex == 1) { # Note that for Cobiac et al parameters, male = 1, female = 0
      sum(take("bmi.sdnl.a"),
          -take("bmi.sdnl.b")*exp(-take("bmi.sdtrend.year.m")*year),
          take("bmi.sdnl.sex")*sex,
          take("bmi.sdtrend.age.m")*age,
          take("bmi.sdnl.age.quad")*age*age,
          take("bmi.sdnl.sexage")*sex*age,
          take("bmi.sdnl.sexagesq")*sex*age*age,
          take("bmi.sdnl.mean")*meanbmi)
    } else if (sex == 0) {
      sum(take("bmi.sdnl.a"),
          -take("bmi.sdnl.b")*exp(-take("bmi.sdtrend.year.f")*year),
          take("bmi.sdnl.sex")*sex,
          take("bmi.sdtrend.age.f")*age,
          take("bmi.sdnl.age.quad")*age*age,
          take("bmi.sdnl.sexage")*sex*age,
          take("bmi.sdnl.sexagesq")*sex*age*age,
          take("bmi.sdnl.mean")*meanbmi)
    } else { print("Error in sex in bmi trend function") }
  } 
} else { print("Error in bestworst")}

# SBP
if (sbpbin==1) { # Temporal trend continues post-2018
  trend.sbp <- function (current.sbp, current.age, sex, current.year) {
    if (current.age<40) {
      new.sbp <- current.sbp
    } else if (sex==1) {
      new.sbp <- sum(current.sbp,
                     take("sbp.trend.age.f"),
                     ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.f"),
                     take("sbp.trend.year.f"),
                     ((current.year-1991+1)^2-(current.year-1991)^2)*take("sbp.trend.year.quad.f"))
    } else if (sex==2) {
      new.sbp <- sum(current.sbp,
                     take("sbp.trend.age.m"),
                     ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.m"),
                     take("sbp.trend.year.m"),
                     ((current.year-1991+1)^2-(current.year-1991)^2)*take("sbp.trend.year.quad.m"))
    } else { print("Error") }
    return(new.sbp)
  }
} else if (sbpbin==0) { # Temporal trend ends in 2018
  trend.sbp <- function (current.sbp, current.age, sex, current.year) {
    if (current.age<40) {
      new.sbp <- current.sbp
    } else if (sex==1) { 
      ifelse(current.year<2019,
             new.sbp <- sum(current.sbp,
                            take("sbp.trend.age.f"),
                            ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.f"),
                            take("sbp.trend.year.f"),
                            ((current.year-1991+1)^2-(current.year-1991)^2)*take("sbp.trend.year.quad.f")),
             new.sbp <- sum(current.sbp,  # Post 2018, no time element in the trend
                            take("sbp.trend.age.f"),
                            ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.f")))
      
    } else if (sex==2) {
      ifelse(current.year<2019,
             new.sbp <- sum(current.sbp,
                            take("sbp.trend.age.m"),
                            ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.m"),
                            take("sbp.trend.year.m"),
                            ((current.year-1991+1)^2-(current.year-1991)^2)*take("sbp.trend.year.quad.m")),
             new.sbp <- sum(current.sbp,  # Post 2018, no time element in the trend
                            take("sbp.trend.age.m"),
                            ((current.age+1)^2-current.age^2)*take("sbp.trend.age.quad.m")))
    } else { print("Error in trend.sbp") }
    return(new.sbp)
  }
} else { print("Error in sbpbin") }

# Cholesterol
if (tchobin==1) { # Temporal trend continues after 2018
  trend.tcho <- function (current.tcho, current.age, sex, current.year) {
    if (sex==1) {
      new.tcho <- sum(current.tcho,
                      take("tcho.trend.age.f"),
                      ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.f"),
                      take("tcho.trend.year.f"),
                      ((current.year-1991+1)^2-(current.year-1991)^2)*take("tcho.trend.year.quad.f"))
    } else if (sex==2) {
      new.tcho <- sum(current.tcho,
                      take("tcho.trend.age.m"),
                      ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.m"),
                      take("tcho.trend.year.m"))
    } else {print("Error")}
    return(new.tcho)
  }
} else if (tchobin==0) { # No temporal trend after 2018
  trend.tcho <- function (current.tcho, current.age, sex, current.year) {
    if (sex==1) {
      ifelse(current.year<2019,
             new.tcho <- sum(current.tcho,
                             take("tcho.trend.age.f"),
                             ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.f"),
                             take("tcho.trend.year.f"),
                             ((current.year-1991+1)^2-(current.year-1991)^2)*take("tcho.trend.year.quad.f")),
             new.tcho <- sum(current.tcho, # Post-2018, no time element in trend
                             take("tcho.trend.age.f"),
                             ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.f")))
    } else if (sex==2) {
      ifelse(current.year<2019,
             new.tcho <- sum(current.tcho,
                             take("tcho.trend.age.m"),
                             ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.m"),
                             take("random.draws$tcho.trend.year.m"),
                             ((current.year-1991+1)^2-(current.year-1991)^2)*take("tcho.trend.year.quad.m")),
             new.tcho <- sum(current.tcho, # Post-2018, no time element in trend
                             take("tcho.trend.age.m"),
                             ((current.age+1)^2-current.age^2)*take("tcho.trend.age.quad.m")))
    } else {print("Error in trend.tcho")}
    return(new.tcho)
  }
} else { print("Error in tchobin") }

## Note that binary risk factors work by first projecting the prevalence for any age and year, 
## and then comparing the chance of risk factor with this prevalence. 
## These functions project the prevalence.
# Smoking
expit = function(z){
  exp(z)/(1+exp(z))
}
if (smokbin==1) { # Temporal trends continue after 2018
  trend.smok.prev <- function (current.age, sex, current.year) {
    if (sex==1) { expit(sum(take("smok.trend.alpha.f"),
                            take("smok.trend.age.f")*current.age,
                            take("smok.trend.age.quad.f")*current.age^2,
                            take("smok.trend.year.f")*(current.year-1991),
                            take("smok.trend.year.quad.f")*(current.year-1991)^2))
    } else if (sex==2) { expit(sum(take("smok.trend.alpha.m"),
                                   take("smok.trend.age.m")*current.age,
                                   take("smok.trend.age.quad.m")*current.age^2,
                                   take("smok.trend.year.m")*(current.year-1991)))
    } else { print("Error") }
  }
} else if (smokbin==0) { # Temporal trends end after 2018
  trend.smok.prev <- function (current.age, sex, current.year) {
    if (sex==1) { ifelse(current.year<2019,
                         expit(sum(take("smok.trend.alpha.f"),
                                   take("smok.trend.age.f")*current.age,
                                   take("smok.trend.age.quad.f")*current.age^2,
                                   take("smok.trend.year.f")*(current.year-1991),
                                   take("smok.trend.year.quad.f")*(current.year-1991)^2)),
                         expit(sum(take("smok.trend.alpha.f"), # Year parameter fixed at 2018
                                   take("smok.trend.age.f")*current.age,
                                   take("smok.trend.age.quad.f")*current.age^2,
                                   take("smok.trend.year.f")*(2018-1991),
                                   take("smok.trend.year.quad.f")*(2018-1991)^2)))
    } else if (sex==2) { ifelse(current.year<2019,
                                expit(sum(take("smok.trend.alpha.m"),
                                          take("smok.trend.age.m")*current.age,
                                          take("smok.trend.age.quad.m")*current.age^2,
                                          take("smok.trend.year.m")*(current.year-1991),
                                          take("smok.trend.year.quad.m")*(current.year-1991)^2)),
                                expit(sum(take("smok.trend.alpha.m"), # Year parameter fixed at 2018
                                          take("smok.trend.age.m")*current.age,
                                          take("smok.trend.age.quad.m")*current.age^2,
                                          take("smok.trend.year.m")*(2018-1991),
                                          take("smok.trend.year.quad.m")*(2018-1991)^2)))
    } else { print("Error in trend.smok.prev") }
  }
} else { print("Error in smokbin") }

# Diabetes
if (diabbin==1) { # Temporal trends continue after 2018
  trend.diab.prev <- function (current.age, sex, current.year) {
    if (sex==1) { expit(sum(take("diab.trend.alpha.f"),
                            take("diab.trend.age.f")*current.age,
                            take("diab.trend.age.quad.f")*current.age^2,
                            take("diab.trend.year.f")*(current.year-1991)))
    } else if (sex==2) { expit(sum(take("diab.trend.alpha.m"),
                                   take("diab.trend.age.m")*current.age,
                                   take("diab.trend.age.quad.m")*current.age^2,
                                   take("diab.trend.year.m")*(current.year-1991),
                                   take("diab.trend.year.quad.m")*(current.year-1991)^2))
    } else { print("Error") }
  }
} else if (diabbin==0) { # Temporal trends end at 2018
  trend.diab.prev <- function (current.age, sex, current.year) {
    if (sex==1) { ifelse(current.year<2019,
                         expit(sum(take("diab.trend.alpha.f"),
                                   take("diab.trend.age.f")*current.age,
                                   take("diab.trend.age.quad.f")*current.age^2,
                                   take("diab.trend.year.f")*(current.year-1991))),
                         expit(sum(take("diab.trend.alpha.f"), # Temporal trends end at 2018
                                   take("diab.trend.age.f")*current.age,
                                   take("diab.trend.age.quad.f")*current.age^2,
                                   take("diab.trend.year.f")*(2018-1991))))
    } else if (sex==2) { ifelse(current.year<2019,
                                expit(sum(take("diab.trend.alpha.m"),
                                          take("diab.trend.age.m")*current.age,
                                          take("diab.trend.age.quad.m")*current.age^2,
                                          take("diab.trend.year.m")*(current.year-1991))),
                                expit(sum(take("diab.trend.alpha.m"), # Temporal trends end at 2018
                                          take("diab.trend.age.m")*current.age,
                                          take("diab.trend.age.quad.m")*current.age^2,
                                          take("diab.trend.year.m")*(2018-1991))))
    } else { print("Error in trend.diab.prev") }
  }
} else { print("Error in diabbin") }

# Function to convert probability to binary estimate for smoking and diabetes
binaryconvert <- function(var1,var2) {
  ifelse(var1>(1-var2),1,0)
}

##############################################################################
#                                                                            #
# 2c. State transition functions                                             #
#                                                                            #
##############################################################################

### p1 transition - probability of new heart attack
inc.fun <- function(a,b,c) {
  min(1,a*(1-b)^(c-2007)/100000)
}

if (treatbin==1) {
  markov.p1 <- function(age,sex,year){ 
    if      (age<18)  { prob <- 0 } # No heart attacks pre-18
    else if (age>100) { prob <- 0 } # No heart attacks post-100 (model kills everyone at 100 anyway) 
    else if (sex==2)  { prob <- inc.fun(a = incdb$baserate.m[age-17],           # age-17 since incdb starts at 18
                                        b = take("treatment.inc.mi.m"),         # Treatment effect centred at 2007
                                        c = year) }                             # (the calibration year)
    else if (sex==1)  { prob <- inc.fun(a = incdb$baserate.f[age-17],
                                        b = take("treatment.inc.mi.f"),
                                        c = year) }
    else              { prob <- c("Error") }
    return(prob)
  }
} else if (treatbin==0) { # This turns off treatment effect for projections past 2018
  markov.p1 <- function(age,sex,year){ 
    if      (age<18)  { prob <- 0 } # No heart attacks pre-18
    else if (age>100) { prob <- 0 } # No heart attacks post-100 (model kills everyone at 100 anyway) 
    else if (sex==2)  { ifelse(year<2019, 
                               prob <- inc.fun(a = incdb$baserate.m[age-17],           # age-17 since incdb starts at 18
                                               b = take("treatment.inc.mi.m"),         # Treatment effect centred at 2007
                                               c = year),
                               prob <- inc.fun(a = incdb$baserate.m[age-17],           # age-17 since incdb starts at 18
                                               b = take("treatment.inc.mi.m"),         # Treatment effect centred at 2007
                                               c = 2018))}                             # (the calibration year)
    else if (sex==1)  { ifelse(year<2019, 
                               prob <- inc.fun(a = incdb$baserate.f[age-17],           # age-17 since incdb starts at 18
                                               b = take("treatment.inc.mi.f"),         # Treatment effect centred at 2007
                                               c = year),
                               prob <- inc.fun(a = incdb$baserate.f[age-17],           # age-17 since incdb starts at 18
                                               b = take("treatment.inc.mi.f"),         # Treatment effect centred at 2007
                                               c = 2018)) }
    else              { prob <- c("Error") }
    return(prob)
  }
} else { print("Error in treatbin") }

### p3: probability of dying from other cause by age and sex
# These parameters are derived in the 'demo.R' script. 

markov.p3 <- function(age,sex,year){
  if      (age<18)  { prob <- 0 } # Under 18s don't die in microPRIME
  else if (sex==2)  { prob <- min(1,exp(take("alpha.mort.m")+take("beta1.mort.m")*age+take("beta2.mort.m")*year)) }
  else if (sex==1)  { prob <- min(1,exp(take("alpha.mort.f")+take("beta1.mort.f")*age+take("beta2.mort.f")*year)) }
  else              { prob <- c("Error") }
  return(prob)
}

### p2: Probability of dying from heart attack by sex, age at heart attack and time since heart attack
# See microPRIME\Data repositories\Parameter space\Functions for case fatality.xlsx for details
# First, calculate probability of 30d case fatality
cf.fun <- function(a,b,c,d,e) {
  min(1,exp(a+b*c)*(1-d)^(e-1999)) # Note, 1999 is the baseline year for the case-fatality functions
}

if(treatbin==1) {
  thirtyd.cf <- function(age,sex,year){ 
    if      (sex==2)  { prob <- cf.fun(a = take("alpha.30dcf.m"),
                                       b = take("beta.30dcf.m"),
                                       c = age,
                                       d = take("treatment.cf.mi.m"),
                                       e = year) } 
    else if (sex==1)  { prob <- cf.fun(a = take("alpha.30dcf.f"),
                                       b = take("beta.30dcf.f"),
                                       c = age,
                                       d = take("treatment.cf.mi.f"),
                                       e = year) }
    else              { prob <- c("Error") }
    return(prob)
  }
} else if(treatbin==0) { # This turns off treatment effect for projections past 2018
  thirtyd.cf <- function(age,sex,year){ 
    if      (sex==2)  { ifelse(year<2019,
                               prob <- cf.fun(a = take("alpha.30dcf.m"),
                                              b = take("beta.30dcf.m"),
                                              c = age,
                                              d = take("treatment.cf.mi.m"),
                                              e = year),
                               prob <- cf.fun(a = take("alpha.30dcf.m"),
                                              b = take("beta.30dcf.m"),
                                              c = age,
                                              d = take("treatment.cf.mi.m"),
                                              e = 2018))
    } 
    else if (sex==1)  { ifelse(year<2019,
                               prob <- cf.fun(a = take("alpha.30dcf.f"),
                                              b = take("beta.30dcf.f"),
                                              c = age,
                                              d = take("treatment.cf.mi.f"),
                                              e = year),
                               prob <- cf.fun(a = take("alpha.30dcf.f"),
                                              b = take("beta.30dcf.f"),
                                              c = age,
                                              d = take("treatment.cf.mi.f"),
                                              e = 2018)) }
    else              { prob <- c("Error") }
    return(prob)
  }
} else { print("Error in treatbin") }

# Now combine 30 day case fatality with background mortality
markov.p2 <- function(age,sex,time.mi,year){ 
  if (time.mi==0) { prob <- min(1,sum(thirtyd.cf(age,sex,year),                                # 30d prob
                                      (1-thirtyd.cf(age,sex,year))*markov.p3(age,sex,year))) } # 30d+ prob
  else            { prob <- markov.p3(age,sex,year) }
  return(prob)
}

