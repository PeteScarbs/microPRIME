### microsim_multiple_popstill_parallel ###

# Last updated: 4th May 2021

# This R script runs the microsimulation model over multiple iterations. Note that a major rewrite took place in
# March 2020, primarily for tidying. A previous version of the R script is available in the SUNDRY folder.

# NOTES: 
# 1. This version of microsim does not allow the population parameters to vary.
# 2. Throughout, female=1, male=2


##############################################################################
#                                                                            #
# Initial settings                                                           #
#                                                                            #
##############################################################################

# Activate packages
#install.packages("foreach")
#install.packages("doParallel")
library("dplyr")
library("splines")
library("foreach")
library("doParallel")

# Model run variables
work        <- "desktop" # Choice between 'desktop' and 'laptop'. Ensures script uses correct directories.
m           <- 450       # The number of iterations of the microsim. This must be less than or equal to the
                         # number of vectors in the 'random draws' file.
n           <- 114000     # 78889 / 114000 The number of agents in each iteration.
y           <- 1998      # The initial year for microsimulation.
l           <- 37        # 20 / 37 The number of years run for each agent in each microsimulation.
core.number <- 18        # The number of cores to run the parallel process on (should be half number available)

# Binary variables to switch on or off effects in projections
bestworst   <- "worst"   # Either "best" or "worst", which use different BMI trend functions that diverge after 2018
treatbin    <- 0         # {0 = No treatment effect after 2018; 1 = Treatment effect after 2018}
sbpbin      <- 0         # {0 = No SBP trend after 2018; 1 = SBP trend after 2018}
tchobin     <- 0         # {0 = No cholesterol trend after 2018; 1 = cholesterol trend after 2018}
smokbin     <- 0         # {0 = No smoking trend after 2018; 1 = smoking trend after 2018}
diabbin     <- 1         # {0 = No diabetes trend after 2018; 1 = Diabetes trend after 2018}

# Input files
if (bestworst == "best") {                   # For varying model parameters.
  random.draw.filename   <- "Random draws_16MAR2021.csv"
} else if (bestworst == "worst") {
  random.draw.filename   <- "Random draws_12APR2021.csv"
  } else { print("Error in random draw filename")}
population.filename    <- "input1_best estimates_13MAR2020_18plus.csv" # Loads the cohort
baserate.filename      <- "input1_best estimatesHSE2006_25FEB2020.csv" # Sets baseline median RR for calibration year.
mi.inc.filename        <- "Smolina_incidence_2007.csv"  # Sets baseline incidence of MI

# Identify final run
final.run   <- "yes"   # If "no" then output years are selected for comparison with external data
                      # and variance is collected. If "yes" then output years for run to final date
                      # are selected and variance is not collected.

# Set age groups for collection of outcome data
a <- c(55, 64)
b <- c(65, 74)
c <- c(75, 84)
d <- c(85, 94)
set <- as.data.frame(rbind(a,b,c,d))

# Set years to collect outcome data
if (final.run == "no") {
  incyears  <- c(1999, 2003, 2007, 2011, 2018) # For emulator output
  prevyears <- c(1999, 2003, 2005, 2006, 2011, 2018)
  eveyears  <- c(1999, 2003, 2007, 2011, 2018)
} else {
  if (final.run == "yes") {
  incyears  <- c(1999, 2003, 2007, 2011, 2015, 2019, 2023, 2027, 2031, 2035) # For final model drawing of results
  prevyears <- c(1999, 2003, 2007, 2011, 2015, 2019, 2023, 2027, 2031, 2035)
  eveyears  <- c(1999, 2003, 2007, 2011, 2015, 2019, 2023, 2027, 2031, 2035)
  } else {"Error"}
}

##############################################################################
#                                                                            #
# Open datasets                                                              #
#                                                                            #
##############################################################################

if (work=="desktop") { 
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Random Draws")
} else {
    if (work=="laptop"){
       setwd("F:\\microPRIME\\Data repositories\\Random Draws")
    }  else { "Error" }
  }
random.draws <- read.csv(random.draw.filename)

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Microsim inputs")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Microsim inputs")
  }  else { "Error" }
}
basedb <- read.csv(baserate.filename)

if (work=="desktop") {
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if(work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  } else { "Error" }
}
mi.db <- read.csv(mi.inc.filename)

if (work=="desktop") {
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Microsim inputs")
} else {
  if(work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Microsim inputs")
  } else { "Error" }
}

pop <- read.csv(population.filename)

##############################################################################
#                                                                            #
# Set baseline MI rates that are constant in all iterations                  #
#                                                                            #
##############################################################################

# Run spline models
spline.ratef <- lm(rate.f ~ bs(age,knots = c(35,54,73)), data = mi.db)
spline.ratem <- lm(rate.m ~ bs(age,knots = c(35,54,73)), data = mi.db)

# Apply to all ages 18 to 100
incdb <- as.data.frame(18:100)
colnames(incdb) <- c("age")
incdb$rate.m <- predict(spline.ratem, newdata = incdb)
incdb$rate.f <- predict(spline.ratef, newdata = incdb)

##############################################################################
#                                                                            #
# Set all functions                                                          #
#                                                                            #
##############################################################################

# There are three sets of function that dictate:
# a. how relative risk of MI for each risk factor varies by sex and age
# b. trends in risk factors over time
# c. state transition between healthy, MI, and death

# All of the functions are in a separate R script (microPRIME functions).

source("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\microsim\\microPRIME functions.R")


##############################################################################
#                                                                            #
# Start the loop around random draw iterations                               #
#                                                                            #
##############################################################################

# Set the parallel backend
c1 <- makeCluster(core.number)  # This nominates a cluster of separate R processes on the different cores
registerDoParallel(c1)          # This starts the parallel process.

model.run.start <- Sys.time()  # Starts the timer for the process

outcomes <- foreach(x = 1:m, 
        .combine = "rbind",
        .packages = "dplyr") %dopar% { 

  ##############################################################################
  #                                                                            #
  # Load all model parameters and decide which should vary between iterations  #
  #                                                                            #
  ##############################################################################
  
  # Set up a dataframe to collect all the parameter information for all of the functions
  df <- as.data.frame(rbind(
    c("minrisk.bmi",         # Parameter name  
      0,                     # 0 = point estimate; 1 = varying parameter
      minrisk.bmi.pe,        # Point estimate value
      999),                  # Varying value (999 if no random draw generated)
    c("sbplow",                   0, sbplow.pe,                999),
    c("tcholow",                  0, tcholow.pe,               999),
    c("combrr.max",               0, combrr.max.pe,            999),
    c("alpha.30dcf.m",            1, alpha.30dcf.m.pe,         random.draws$alpha.30dcf.m[x]),
    c("beta.30dcf.m",             0, beta.30dcf.m.pe,          999), 
    c("alpha.30dcf.f",            1, alpha.30dcf.f.pe,         random.draws$alpha.30dcf.f[x]),
    c("beta.30dcf.f",             0, beta.30dcf.f.pe,          999), 
    c("alpha.mort.m",             0, alpha.mort.m.pe,          999),
    c("beta1.mort.m",             0, beta1.mort.m.pe,          999),
    c("beta2.mort.m",             0, beta2.mort.m.pe,          999),
    c("alpha.mort.f",             0, alpha.mort.f.pe,          999),
    c("beta1.mort.f",             0, beta1.mort.f.pe,          999),
    c("beta2.mort.f",             0, beta2.mort.f.pe,          999),
    c("treatment.inc.mi.m",       0, treatment.inc.mi.m.pe,    999),
    c("treatment.cf.mi.m",        1, treatment.cf.mi.m.pe,     random.draws$treatment.cf.mi.m[x]),
    c("treatment.inc.mi.f",       0, treatment.inc.mi.f.pe,    999), 
    c("treatment.cf.mi.f",        1, treatment.cf.mi.f.pe,     random.draws$treatment.cf.mi.f[x]),
    c("bmi.m.r1525",              1, bmi.m.r1525.pe,           random.draws$bmi.m.r1525[x]),
    c("bmi.m.r2550",              1, bmi.m.r2550.pe,           random.draws$bmi.m.r2550[x]),
    c("bmi.f.r1525",              1, bmi.f.r1525.pe,           random.draws$bmi.f.r1525[x]),
    c("bmi.f.r2550",              1, bmi.f.r2550.pe,           random.draws$bmi.f.r2550[x]),
    c("sbp.alpha.m",              1, sbp.alpha.m.pe,           random.draws$sbp.alpha.m[x]),
    c("sbp.beta.m",               0, sbp.beta.m.pe,            999), 
    c("sbp.alpha.f",              1, sbp.alpha.f.pe,           random.draws$sbp.alpha.f[x]),
    c("sbp.beta.f",               0, sbp.beta.f.pe,            999),
    c("tcho.alpha.m",             1, tcho.alpha.m.pe,          random.draws$tcho.alpha.m[x]),
    c("tcho.beta.m",              0, tcho.beta.m.pe,           999), 
    c("tcho.alpha.f",             1, tcho.alpha.f.pe,          random.draws$tcho.alpha.f[x]),
    c("tcho.beta.f",              0, tcho.beta.f.pe,           999),
    c("smok.alpha.m",             1, smok.alpha.m.pe,          random.draws$smok.alpha.m[x]),
    c("smok.beta.m",              0, smok.beta.m.pe,           999), 
    c("smok.alpha.f",             1, smok.alpha.f.pe,          random.draws$smok.alpha.f[x]),
    c("smok.beta.f",              0, smok.beta.f.pe,           999),
    c("diab.m",                   1, diab.m.pe,                random.draws$diab.m[x]),
    c("diab.f",                   1, diab.f.pe,                random.draws$diab.f[x]),
    c("prev.alpha.m",             1, prev.alpha.m.pe,          random.draws$prev.alpha.m[x]), 
    c("prev.alpha.f",             1, prev.alpha.f.pe,          random.draws$prev.alpha.f[x]),
    c("prev.beta.age.m",          0, prev.beta.age.m.pe,       999),
    c("prev.beta.age.f",          0, prev.beta.age.f.pe,       999),
    c("prev.beta.time.m",         0, prev.beta.time.m.pe,      999),
    c("prev.beta.time.f",         0, prev.beta.time.f.pe,      999),
    c("sbp.trend.age.m",          1, sbp.trend.age.m.pe,       random.draws$sbp.trend.age.m[x]),
    c("sbp.trend.age.f",          1, sbp.trend.age.f.pe,       random.draws$sbp.trend.age.f[x]),
    c("sbp.trend.age.quad.m",     0, sbp.trend.age.quad.m.pe,  999),
    c("sbp.trend.age.quad.f",     0, sbp.trend.age.quad.f.pe,  999),
    c("sbp.trend.year.m",         1, sbp.trend.year.m.pe,      random.draws$sbp.trend.year.m[x]),
    c("sbp.trend.year.f",         1, sbp.trend.year.f.pe,      random.draws$sbp.trend.year.f[x]),
    c("sbp.trend.year.quad.m",    0, sbp.trend.year.quad.m.pe, 999),
    c("sbp.trend.year.quad.f",    0, sbp.trend.year.quad.f.pe, 999),
    c("tcho.trend.age.m",         1, tcho.trend.age.m.pe,      random.draws$tcho.trend.age.m[x]),
    c("tcho.trend.age.f",         1, tcho.trend.age.f.pe,      random.draws$tcho.trend.age.f[x]),
    c("tcho.trend.age.quad.m",    0, tcho.trend.age.quad.m.pe, 999),
    c("tcho.trend.age.quad.f",    0, tcho.trend.age.quad.f.pe, 999),
    c("tcho.trend.year.m",        1, tcho.trend.year.m.pe,     random.draws$tcho.trend.year.m[x]),
    c("tcho.trend.year.f",        1, tcho.trend.year.f.pe,     random.draws$tcho.trend.year.f[x]),
    c("tcho.trend.year.quad.f",   0, tcho.trend.year.quad.f.pe,999),
    c("smok.trend.alpha.m",       0, smok.trend.alpha.m.pe,    999), 
    c("smok.trend.age.m",         1, smok.trend.age.m.pe,      random.draws$smok.trend.age.m[x]),
    c("smok.trend.age.quad.m",    0, smok.trend.age.quad.m.pe, 999),
    c("smok.trend.year.m",        1, smok.trend.year.m.pe,     random.draws$smok.trend.year.m[x]),
    c("smok.trend.alpha.f",       0, smok.trend.alpha.f.pe,    999),
    c("smok.trend.age.f",         1, smok.trend.age.f.pe,      random.draws$smok.trend.age.f[x]),
    c("smok.trend.age.quad.f",    0, smok.trend.age.quad.f.pe, 999),
    c("smok.trend.year.f",        1, smok.trend.year.f.pe,     random.draws$smok.trend.year.f[x]),
    c("smok.trend.year.quad.f",   0, smok.trend.year.quad.f.pe,999),
    c("diab.trend.alpha.m",       0, diab.trend.alpha.m.pe,    999),
    c("diab.trend.age.m",         1, diab.trend.age.m.pe,      random.draws$diab.trend.age.m[x]),
    c("diab.trend.age.quad.m",    0, diab.trend.age.quad.m.pe, 999),
    c("diab.trend.alpha.f",       0, diab.trend.alpha.f.pe,    999),
    c("diab.trend.age.f",         1, diab.trend.age.f.pe,      random.draws$diab.trend.age.f[x]),
    c("diab.trend.age.quad.f",    0, diab.trend.age.quad.f.pe, 999),
    c("diab.trend.year.m",        1, diab.trend.year.m.pe,     random.draws$diab.trend.year.m[x]),
    c("diab.trend.year.quad.m",   0,diab.trend.year.quad.m.pe, 999),
    c("diab.trend.year.f",        1, diab.trend.year.f.pe,     random.draws$diab.trend.year.f[x]),
    c("bmi.meantrend.const",      0, bmi.meantrend.const.pe,   999),
    c("bmi.meantrend.sex",        0, bmi.meantrend.sex.pe,     999),
    c("bmi.trend.age.m",          1, bmi.meantrend.age.pe,     random.draws$bmi.trend.age.m[x]),
    c("bmi.trend.age.f",          1, bmi.meantrend.age.pe,     random.draws$bmi.trend.age.f[x]),
    c("bmi.trend.year.m",         1, bmi.meantrend.year.pe,    random.draws$bmi.trend.year.m[x]),
    c("bmi.trend.year.f",         1, bmi.meantrend.year.pe,    random.draws$bmi.trend.year.f[x]),
    c("bmi.meantrend.age.quad",   0,bmi.meantrend.age.quad.pe, 999),
    c("bmi.meantrend.year.quad",  0,bmi.meantrend.year.quad.pe,999),  
    c("bmi.sdtrend.const",        0, bmi.sdtrend.const.pe,     999),
    c("bmi.sdtrend.sex",          0, bmi.sdtrend.sex.pe,       999),
    c("bmi.sdtrend.mean",         0, bmi.sdtrend.mean.pe,      999),
    c("bmi.sdtrend.mean.quad",    0, bmi.sdtrend.mean.quad.pe, 999),
    c("bmi.sdtrend.age.m",        0, bmi.sdtrend.age.pe,       999),
    c("bmi.sdtrend.age.f",        0, bmi.sdtrend.age.pe,       999),
    c("bmi.sdtrend.year.m",       0, bmi.sdtrend.year.pe,      999),
    c("bmi.sdtrend.year.f",       0, bmi.sdtrend.year.pe,      999),
    c("bmi.sdtrend.age.quad",     0, bmi.sdtrend.age.quad.pe,  999),
    c("bmi.sdtrend.year.quad",    0, bmi.sdtrend.year.quad.pe, 999),
    c("bmi.sdtrend.sexage",       0, bmi.sdtrend.sexage.pe,    999),
    c("bmi.sdtrend.sexagesq",     0, bmi.sdtrend.sexagesq.pe,  999),
    c("bmi.sdtrend.ageyear",      0, bmi.sdtrend.ageyear.pe,   999),
    c("bmi.sdtrend.ageyearsq",    0, bmi.sdtrend.ageyearsq.pe, 999),
    c("bmi.meannl.a",             0, bmi.meannl.a.pe,          999),
    c("bmi.meannl.b",             0, bmi.meannl.b.pe,          999),
    c("bmi.meannl.year.m",        0, bmi.meannl.year.pe,       999),
    c("bmi.meannl.year.f",        0, bmi.meannl.year.pe,       999),
    c("bmi.meannl.sex",           0, bmi.meannl.sex.pe,        999),
    c("bmi.meannl.age.m",         0, bmi.meannl.age.pe,        999),
    c("bmi.meannl.age.f",         0, bmi.meannl.age.pe,        999),
    c("bmi.meannl.age.quad",      0, bmi.meannl.age.quad.pe,   999),
    c("bmi.meannl.sexage",        0, bmi.meannl.sexage.pe,     999),
    c("bmi.meannl.sexagesq",      0, bmi.meannl.sexagesq.pe,   999),
    c("bmi.sdnl.a",               0, bmi.sdnl.a.pe,            999),
    c("bmi.sdnl.b",               0, bmi.sdnl.b.pe,            999),
    c("bmi.sdnl.year.m",          0, bmi.sdnl.year.pe,         999),
    c("bmi.sdnl.year.f",          0, bmi.sdnl.year.pe,         999),
    c("bmi.sdnl.sex",             0, bmi.sdnl.sex.pe,          999),
    c("bmi.sdnl.age.m",           0, bmi.sdnl.age.pe,          999),
    c("bmi.sdnl.age.f",           0, bmi.sdnl.age.pe,          999),
    c("bmi.sdnl.age.quad",        0, bmi.sdnl.age.quad.pe,     999),
    c("bmi.sdnl.sexage",          0, bmi.sdnl.sexage.pe,       999),
    c("bmi.sdnl.sexagesq",        0, bmi.sdnl.sexagesq.pe,     999),
    c("bmi.sdnl.mean",            0, bmi.sdnl.mean.pe,         999)))
  
  colnames(df) <- c("name","onoff","pointestimate","varyingparameter")
  df.name <- as.character(df$name)
  df$onoff <- as.integer(as.character(df$onoff))
  df$pointestimate <- as.numeric(as.character(df$pointestimate))
  df$varyingparameter <- as.numeric(as.character(df$varyingparameter))
  
  # Add a column to draw the final model parameter from
  decide <- function(p,q,r) {
    q + p*(r - q)
  }
  df$modelparameter <- mapply(p = df$onoff,
                              q = df$pointestimate,
                              r = df$varyingparameter,
                              FUN = decide)
                  
  ##############################################################################
  #                                                                            #
  # Set up baseline dataset for calibrating incidence rates                    #
  #                                                                            #
  ##############################################################################
  
  # Convert smoking and diabetes from probability to binary variable
  basedb$smok.prev <- mapply(current.age = basedb$age, 
                             sex = basedb$sex2, 
                             current.year = 2006,   # Same as baseline rate cohort.
                             FUN = trend.smok.prev)
  basedb$diab.prev <- mapply(current.age = basedb$age, 
                             sex = basedb$sex2, 
                             current.year = 2006,
                             FUN = trend.diab.prev)
  
  basedb$smok <- mapply(binaryconvert, basedb$smok.p, basedb$smok.prev)
  basedb$diab <- mapply(binaryconvert, basedb$diab.p, basedb$diab.prev)
  
    basedb$rr <- mapply(sex = basedb$sex2,   # This calculates the combined RR for the basline rate cohort
                      age = basedb$age,
                      bmi = basedb$bmi,
                      sbp = basedb$sbp,
                      tcho = basedb$tcho,
                      smok = basedb$smok,
                      diab = basedb$diab,
                      FUN = rr.fun)
  
  # Collapse on median rr levels by age category and sex
  baseage <- basedb %>% group_by(age, sex2) %>% summarise(rr = median(rr))  
  
  # Generate data frame with baseline rates for age 18 to 100 both sexes
  baserate.function <- function (unadj.rate, rr.adj) {
    max(0,unadj.rate) / max(1,rr.adj)
  }
  
  incdb$baserate.m <- mapply(unadj.rate = incdb$rate.m,
                             rr.adj = baseage$rr[baseage$sex2==2],
                             FUN = baserate.function)
  incdb$baserate.f <- mapply(unadj.rate = incdb$rate.f,
                             rr.adj = baseage$rr[baseage$sex2==1],
                             FUN = baserate.function)

  ##############################################################################
  #                                                                            #
  # Iterate around agents in cohort                                            #
  #                                                                            #
  ##############################################################################

  ## Set the empty dataframe for collecting  the microsim output
  # This is what will be built, row by row, by iterating around agents.
  # Each agent generates a temporary data.frame (see below), which is reduced to a single row of data for
  # collection of outcomes.
  micro.db <- data.frame(id=integer(),        # Agent id number
                         sex=integer(),       # {1 = Female; 2 = Male}
                         age=integer(),       # Age at initial year
                         bmi=numeric(),       # BMI (kg/m2) at initial year
                         sbp=numeric(),       # SBP (mmHg) at initial year
                         tcho=numeric(),      # Total cholesterol (mmol/l) at initial year
                         smok=integer(),      # {1 = smoker at initial year; 0 = Not}
                         diab=integer(),      # {1 = Has diabetes at initial year; 0 = Not}
                         year.mi1=integer(),  # {0 = No first MI; other = year of first MI}
                         year.mi2=integer(),  # {0 = No second MI; other = year of second MI}
                         year.mi3=integer(),  # {0 = No third MI; other = year of third MI}
                         year.mi4=integer(),  # {0 = No fourth MI; other = year of foyrth MI}
                         year.dead=integer()) # {0 = No death; other = year of death}
  
  # Begin loop around agents
  for (j in 1:n){
  
    # First for each agent in pop, set an initial dataset of relevant characteristics
    y0              <- c(y, pop[j,3:10])
    micro           <- data.frame(y0)
    colheader       <- c("year", "sex", "age", "bmi", "sbp", "tcho", "smok.p", "diab.p", "prev.mi")
    colnames(micro) <- colheader
    
    # Then run this dataset of characteristics forward for the number of years specified in the model.
    for (i in 1:l) {
      newyear <- c(micro$year[i]+1,
                   micro$sex[1],
                   micro$age[i]+1,
                   micro$bmi[1], # Holding for now
                   trend.sbp (micro$sbp[i], micro$age[i],micro$sex[1],micro$year[i]),
                   trend.tcho(micro$tcho[i],micro$age[i],micro$sex[1],micro$year[i]),
                   micro$smok.p[1],
                   micro$diab.p[1],
                   micro$prev.mi[1])
      micro <- rbind(micro,newyear)
    } 
    
    # Add columns for prevalence of smoking and diabetes (for conversion of smok.p and diab.p to binary)
    micro$smok.prev <- mapply(trend.smok.prev, micro$age, micro$sex, micro$year)
    micro$diab.prev <- mapply(trend.diab.prev, micro$age, micro$sex, micro$year)
    
    # Convert probability of smoking and diabetes into binary estimates
    micro$smok <- mapply(binaryconvert, micro$smok.p, micro$smok.prev)
    micro$diab <- mapply(binaryconvert, micro$diab.p, micro$diab.prev)
    
    # Add column for projected mean and sd of lnBMI
    micro$mean.lnbmi <- mapply(age = micro$age,
                               sex = micro$sex - 1, #NB: adjustment to Cobiac et al. parameters
                               year = micro$year - 2003,
                               FUN = mean.bmi.trend)
    micro$sd.lnbmi <- mapply(age = micro$age,
                             sex = micro$sex - 1, 
                             year = micro$year - 2003,
                             meanbmi = micro$mean.lnbmi,
                             FUN = sd.bmi.trend)
    
    # Add column for lnbmi z-score
    micro$bmi.z <- z.fun(bmi = micro$bmi[1],
                         meanlnbmi = micro$mean.lnbmi[1],
                         sdlnbmi = micro$sd.lnbmi[1])
    
    # Convert BMI to projection
    micro$bmi <- mapply(z = micro$bmi.z,
                        mean = micro$mean.lnbmi,
                        sd = micro$sd.lnbmi,
                        FUN = backz)
    
    ### Add incident MIs
    ## Generate two column vectors, that will be added to micro once filled:
    ## mi {0 = no MI that year; 1 = MI that year}
    ## time.mi, continuous variable counting years since last MI
    # First year of column vectors- assume that those who have had previous MI that it happened 5 years before
    # start of model
    if (micro$prev.mi[1]==1) { 
      mi      <- 0             
      time.mi <- 5             
    } else {
      mi      <- 0             
      time.mi <- 0
    }
    
    # For remaining years, incident MIs are based on RRs and p1 
    for (i in 2:(l+1)) {
      newmi <- rbinom(1,1,min(1, # 'min' needed so probability is between 0 and 1
                              min(take("combrr.max"), rrbmi(micro$sex[i],micro$bmi[i])*   
                                              rrsbp(micro$age[i],micro$sex[i],micro$sbp[i])*   
                                              rrtcho(micro$age[i],micro$sex[i],micro$tcho[i])* 
                                              rrsmok(micro$age[i],micro$sex[i],micro$smok[i])*
                                              rrdiab(micro$sex[i],micro$diab[i])*
                                              rrprev(micro$age[i],micro$sex[i],time.mi[i-1]))*
                              markov.p1(micro$age[i],micro$sex[i],micro$year[i])))
      if (newmi==1) {
        newtime.mi <- 1
      } else if (time.mi[i-1]>0) {
        newtime.mi <- time.mi[i-1]+1
      } else {
        newtime.mi <- 0
      }
      mi      <- rbind(mi, newmi)
      time.mi <- rbind(time.mi, newtime.mi)
      row.names(mi)[i]      <- i 
      row.names(time.mi)[i] <- i
    }
    
    # Combine with micro
    mi                <- data.frame(mi)
    time.mi           <- data.frame(time.mi)
    colmi             <- c("mi")
    colnames(mi)      <- colmi
    coltime.mi        <- c("time.mi")
    colnames(time.mi) <- coltime.mi
    micro             <- cbind(micro,mi,time.mi)
    
    # Generate new column vector to add to micro:
    # year.mi, the year when the most recent mi occurred.
    # Generate first year
    if (micro$prev.mi[1]==1) { 
      year.mi <- y-5           
    } else {
      year.mi <- 0
    }
    
    # Subsequent years
    for (i in 2:(l+1)) {
      newyear.mi <- if (micro$mi[i]==0) {
        year.mi[i-1]
      } else { 
        micro$mi[i]*micro$year[i]
        }
      year.mi <- rbind(year.mi, newyear.mi)
      row.names(year.mi)[i] <- i
    }
    
    # Combine with micro
    year.mi           <- data.frame(year.mi)
    colyearmi         <- c("year.mi") 
    colnames(year.mi) <- colyearmi
    micro             <- cbind(micro,year.mi)
    
    # Generate new column vector to add to micro:
    # count.mi, a running total of the number of MIs an agent has had.
    # Generate first year
    if (micro$prev.mi[1]==1) { 
      count.mi <- 1            
    } else {
      count.mi <- 0
    }
    
    # Subsequent years
    for (i in 2:(l+1)) {
      newcount.mi <- if (micro$mi[i]==0) {
        count.mi[i-1]
      } else {
        count.mi[i-1]+1
      }
      count.mi <- rbind(count.mi, newcount.mi)
      row.names(count.mi)[i] <- i
    }
    
    # Combine with micro
    count.mi           <- data.frame(count.mi)
    colcountmi         <- c("count.mi") 
    colnames(count.mi) <- colcountmi
    micro              <- cbind(micro,count.mi)
    
    # Generate new column vector to add to micro:
    # death, {0 = alive; other = year of death}
    # Generate first year, where by definition all agents are alive
    death <- 0
    
    # Subsequent years
    for (i in 2:(l+1)) {
      newdeath <- if (death[i-1]>0) { # If someone is dead, they can't die again!
        death[i-1]
        }
      else if (micro$count.mi[i]>=4) { # This ensures that a fourth MI kills the agent
        micro$year[i]
      }
      else if (micro$year.mi[i]>0) { # If agent has had MI, then p2 controls mortality
        micro$year[i]*rbinom(1,1,min(1,markov.p2(age = micro$age[i],
                                                 sex = micro$sex[1],
                                                 time.mi = micro$time.mi[i]-1, # -1 is necessary so that in the year of an MI, time.mi==0 (see microPRIME functions)
                                                 year = micro$year[i])))
        }
      else {                         # If no MI, then p3 controls mortality
        micro$year[i]*rbinom(1,1,min(1,markov.p3(age = micro$age[i],
                                                 sex = micro$sex[i],
                                                 year = micro$year[i])))
        }
      death <- rbind(death, newdeath)
      row.names(death)[i] <- i
    }
    
    # Combine with micro
    death           <- data.frame(death)
    coldeath        <- c("year.death")
    colnames(death) <- coldeath
    micro           <- cbind(micro,death)
    
    # Extracting and storing the important information
    micro.db.j <- c(j,
                  micro$sex[1],
                  micro$age[1],
                  micro$bmi[1],
                  micro$sbp[1],
                  micro$tcho[1],
                  micro$smok[1],
                  micro$diab[1],
                  if      (max(micro$year.death)==0)                 {max(micro$year.mi[micro$count.mi<=1])}
                  else if (max(micro$year.mi[micro$count.mi<=1])>max(micro$year.death))                  {0}
                  else                                               {max(micro$year.mi[micro$count.mi<=1])},
                  if      (max(micro$count.mi)<=1)                                                       {0}
                  else if (max(micro$year.death)==0)                 {max(micro$year.mi[micro$count.mi==2])}
                  else if (max(micro$year.mi[micro$count.mi==2])>max(micro$year.death))                  {0}
                  else                                               {max(micro$year.mi[micro$count.mi==2])},
                  if      (max(micro$count.mi)<=2)                                                       {0}
                  else if (max(micro$year.death)==0)                 {max(micro$year.mi[micro$count.mi==3])}
                  else if (max(micro$year.mi[micro$count.mi==3])>max(micro$year.death))                  {0}
                  else                                               {max(micro$year.mi[micro$count.mi==3])},
                  if      (max(micro$count.mi)<=3)                                                       {0}
                  else if (max(micro$year.death)==0)                 {max(micro$year.mi[micro$count.mi==4])}
                  else if (max(micro$year.mi[micro$count.mi==4])>max(micro$year.death))                  {0}
                  else                                               {max(micro$year.mi[micro$count.mi==4])},
                  max(micro$year.death))
    micro.db.j <- data.frame(rbind(micro.db.j))
    colmicro.db <- c("id", "sex", "age", "bmi", "sbp", "tcho", "smok", "diab", 
                     "year.mi1", "year.mi2", "year.mi3", "year.mi4", "year.dead")
    colnames(micro.db.j) <- colmicro.db
    
    # Store in micro.db
    micro.db <- rbind(micro.db, micro.db.j)
    
  }
  
  ##############################################################################
  #                                                                            #
  # Collect output data                                                        #
  #                                                                            #
  ##############################################################################
  
  # Add a constant for counting agents
  micro.db$cons <- 1
  
  # Differentiate between final run and training data run
  if (final.run == "no") {
    ## Extract incidence rates per 100,000 and variance (NB: first MI only)
    # Loop round age groups in set
    for (j in incyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the incidence outcomes
        inc.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & micro.db$year.mi1==j])
        
        inc.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & micro.db$year.mi1==j])
        
        assign(paste0("incrate.m.",j,".",i),inc.m*100000/num.m)
        assign(paste0("incrate.f.",j,".",i),inc.f*100000/num.f)
        
        # Variance is based on square of the standard error of a propotion
        assign(paste0("var.inc.m.",j,".",i),
               get(paste0("incrate.m.",j,".",i))*(100000-get(paste0("incrate.m.",j,".",i)))/(100000*num.m))  
        assign(paste0("var.inc.f.",j,".",i),
               get(paste0("incrate.f.",j,".",i))*(100000-get(paste0("incrate.f.",j,".",i)))/(100000*num.f))
        
      } 
      
    }  
    
    # Extract prevalence data
    for (j in prevyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the prevalence outcomes
        prev.m <- sum(micro.db$cons[micro.db$sex==2
                                    & micro.db$age+(j-y)>=set[i,1]
                                    & micro.db$age+(j-y)<=set[i,2]
                                    & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                    & micro.db$year.mi1!=0
                                    & micro.db$year.mi1<=j])
        
        prev.f <- sum(micro.db$cons[micro.db$sex==1
                                    & micro.db$age+(j-y)>=set[i,1]
                                    & micro.db$age+(j-y)<=set[i,2]
                                    & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                    & micro.db$year.mi1!=0
                                    & micro.db$year.mi1<=j])
        
        # Extract prevalence rates as proportions
        assign(paste0("prevrate.m.",j,".",i),prev.m/num.m)
        assign(paste0("prevrate.f.",j,".",i),prev.f/num.f)
        
        # Variance is based on square of the standard error of a propotion
        assign(paste0("var.prev.m.",j,".",i),
               get(paste0("prevrate.m.",j,".",i))*(1-get(paste0("prevrate.m.",j,".",i)))/num.m)
        assign(paste0("var.prev.f.",j,".",i),
               get(paste0("prevrate.f.",j,".",i))*(1-get(paste0("prevrate.f.",j,".",i)))/num.f)
        
      }
      
    }
    
    # Extract event rate data (i.e rate of any heart attack, not just first)
    # Loop round age groups in set
    for (j in eveyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the event rate outcomes
        eve.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & (micro.db$year.mi1==j | micro.db$year.mi2==j | micro.db$year.mi3==j | micro.db$year.mi4==j )])
        
        eve.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & (micro.db$year.mi1==j | micro.db$year.mi2==j | micro.db$year.mi3==j | micro.db$year.mi4==j )])
        
        # Extract event rates per 100,000
        assign(paste0("everate.m.",j,".",i),eve.m*100000/num.m) 
        assign(paste0("everate.f.",j,".",i),eve.f*100000/num.f)
        
        # Variance is based on square of standard error of proportion
        assign(paste0("var.eve.m.",j,".",i),
               get(paste0("everate.m.",j,".",i))*(100000-get(paste0("everate.m.",j,".",i)))/(100000*num.m))  
        assign(paste0("var.eve.f.",j,".",i),
               get(paste0("everate.f.",j,".",i))*(100000-get(paste0("everate.f.",j,".",i)))/(100000*num.f))
        
      } 
      
    }  
    
    # Collect data into single dataset
    end <- NULL
    
    for (j in incyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("incrate.m.",j,".",i)),get(paste0("var.inc.m.",j,".",i)))
      }
    }
    
    for (j in incyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("incrate.f.",j,".",i)),get(paste0("var.inc.f.",j,".",i)))
      }
    }
    
    for (j in prevyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("prevrate.m.",j,".",i)),get(paste0("var.prev.m.",j,".",i)))
      }
    }
    
    for (j in prevyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("prevrate.f.",j,".",i)),get(paste0("var.prev.f.",j,".",i)))
      }
    }
    
    for (j in eveyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("everate.m.",j,".",i)),get(paste0("var.eve.m.",j,".",i)))
      }
    }
    
    for (j in eveyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("everate.f.",j,".",i)),get(paste0("var.eve.f.",j,".",i)))
      }
    }
    
    # outcomes <- rbind(outcomes, end)
    
    # Loop finished - collect data
    return(cbind(x,end))
  
  } else { # Code below is what happens if it is the final run
  
    ## Extract incidence rates per 100,000 (NB: first MI only)
    # Loop round age groups in set
    for (j in incyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the incidence outcomes
        inc.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & micro.db$year.mi1==j])
        
        inc.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & micro.db$year.mi1==j])
        
        assign(paste0("incrate.m.",j,".",i),inc.m*100000/num.m)
        assign(paste0("incrate.f.",j,".",i),inc.f*100000/num.f)
      } 
      
    }  
    
    # Extract prevalence data
    for (j in prevyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the prevalence outcomes
        prev.m <- sum(micro.db$cons[micro.db$sex==2
                                    & micro.db$age+(j-y)>=set[i,1]
                                    & micro.db$age+(j-y)<=set[i,2]
                                    & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                    & micro.db$year.mi1!=0
                                    & micro.db$year.mi1<=j])
        
        prev.f <- sum(micro.db$cons[micro.db$sex==1
                                    & micro.db$age+(j-y)>=set[i,1]
                                    & micro.db$age+(j-y)<=set[i,2]
                                    & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                    & micro.db$year.mi1!=0
                                    & micro.db$year.mi1<=j])
        
        # Extract prevalence rates as proportions
        assign(paste0("prevrate.m.",j,".",i),prev.m/num.m)
        assign(paste0("prevrate.f.",j,".",i),prev.f/num.f)
      }
      
    }
    
    # Extract event rate data (i.e rate of any heart attack, not just first)
    # Loop round age groups in set
    for (j in eveyears) {
      for (i in 1:4) {
        # First get the denominator, which is the number of people in age range alive in outcome year
        num.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        num.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j) ])
        
        # Numerator for the event rate outcomes
        eve.m <- sum(micro.db$cons[micro.db$sex==2
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & (micro.db$year.mi1==j | micro.db$year.mi2==j | micro.db$year.mi3==j | micro.db$year.mi4==j )])
        
        eve.f <- sum(micro.db$cons[micro.db$sex==1
                                   & micro.db$age+(j-y)>=set[i,1]
                                   & micro.db$age+(j-y)<=set[i,2]
                                   & (micro.db$year.dead==0 | micro.db$year.dead>=j)
                                   & (micro.db$year.mi1==j | micro.db$year.mi2==j | micro.db$year.mi3==j | micro.db$year.mi4==j )])
        
        # Extract event rates per 100,000
        assign(paste0("everate.m.",j,".",i),eve.m*100000/num.m) 
        assign(paste0("everate.f.",j,".",i),eve.f*100000/num.f)
      } 
      
    }  
    
    # Collect data into single dataset
    end <- NULL
    
    for (j in incyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("incrate.m.",j,".",i)))
      }
    }
    
    for (j in incyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("incrate.f.",j,".",i)))
      }
    }
    
    for (j in prevyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("prevrate.m.",j,".",i)))
      }
    }
    
    for (j in prevyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("prevrate.f.",j,".",i)))
      }
    }
    
    for (j in eveyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("everate.m.",j,".",i)))
      }
    }
    
    for (j in eveyears) {
      for (i in 1:4) {
        end <- cbind(end,get(paste0("everate.f.",j,".",i)))
      }
    }
    
    # outcomes <- rbind(outcomes, end)
    
    # Loop finished - collect data
    return(cbind(x,end))
  }
} # This ends the foreach loop

# Find model run time

stopCluster(c1) # This stops the parallel process.

model.run.end <- Sys.time() 
model.run.end - model.run.start
 

# Save model outcomes

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Training data")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Training data")
  }  else { "Error" }
}

write.csv(outcomes, file = "Model 18 WORST CASE_04MAY2021.csv") # Change this name each time.

