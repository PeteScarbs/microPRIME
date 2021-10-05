### agent fill ###

# Last updated: 13th March 2020

# Note - this script can be used as a standalone to generate a population csv that is then used in microsim. If 
# population parameters are allowed to vary, then most of this script is incorporated into microsim anyway.

# This R script generates an artificial population with a joint distribution of the following
# variables:
# 1. Age (NB: under 18s are new entries that won't be included in the model until after the starting year)
# 2. Sex
# 3. Percentile risk of current smoking
# 4. BMI
# 5. Systolic blood pressure
# 6. Total blood cholesterol
# 7. Percentile risk of diabetes
# 8. Previous MI (0 = no; 1 = yes)
# with joint distribution described by plausible parameters drawn from the 'random draws' data
# repository and the population stratfied by sex.

# Also generates whether or not the agent has migrated to the model yet, and what year they will migrate to the model.

# Load necessary packages
# install.packages("copula")
# install.packages("rms")
# install.packages("data.table")
# install.packages("plyr")
library("copula")
library("rms")
library("data.table") 
library("plyr")

# model.run.start <- Sys.time()

work <- "desktop" # Choice between 'desktop' and 'laptop'. Ensures script uses correct directories.

# Set number of agents 
# (NB: this should be half of n in microsim. Half because it generates n men and n women and combines for full sample)
n <- 57000

# Set number of datasets that need to be produced (NB: this should match m in microsim.R and n in parameter draw.R)
m <- 1 # This should be the number of iterations across populations.

# Set minimum and maximum age for the generated population. For calibration of microsim, this should be 55 - (difference between 2012 and baseline year), 
# as no events will be recorded for people under the age of 55 in the model run.
min.age <- 18
max.age <- 95

# Open the random draws data repository
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}
# NB: change this if drawing from 'Random draws'

vector.sample <- read.csv("Parameter space_best estimates_28NOV2019.csv") # NB: To update for each run


## Constructing the age structure of the population
# Step 1: Intro the population of the baseline year

yr <- 1998 # Set the baseline year. Should be the same as the baseline year in risk factors.R

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Demographic trends")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Demographic trends")
  }  else { "Error" }
}

population <- read.csv("pop.mort04AUG2017.csv") # To update the input file, need to run 'demo_preparing dataset.R'

# Reduce to baseline population data only
pop.base <- population[population$year==yr,2:5]

# Step 2: Divide the 85+ population into single age categories
# Plot population against age by sex
par(mfrow=c(2,1))
plot(pop.base$age[pop.base$sex==1],pop.base$population[pop.base$sex==1],
     xlab = "Age", ylab = "Population, men")
plot(pop.base$age[pop.base$sex==2],pop.base$population[pop.base$sex==2],
     xlab = "Age", ylab = "Population, women")
# Reduction in population by age between late 70s and mid 80s looks linear. 

# Produce a linear trend with the following constraints:
# 1. 85+ group is split between 85 to 100 year olds
# 2. Total sum of population in 85 to 100 groups equals 85+ figure
# 3. Amount for 85 year olds is a linear continuation of trend from 80-84 year olds

# A bit of maths, shows that this trend can be estimated as y = mx + c, where:
# x = 0, when age is 85
# c = extrapolated trend from 80-84 year olds
# m = (X - 16c)/120, where X = population in 85+ group
# NB: This works - but can produce negative numbers in some groups. So run to 100, and then cut back to sensible
# number to avoid negatives.

# Calculate c for men and women
extrap.func <- function (a, b, num.between, extrap.num) {
  b + ((b-a)/num.between)*extrap.num
}
c.m <- extrap.func(a = pop.base$population[pop.base$sex==1&pop.base$age==80],
                   b = pop.base$population[pop.base$sex==1&pop.base$age==84],
                   num.between = 4,
                   extrap.num = 1)
c.f <- extrap.func(a = pop.base$population[pop.base$sex==2&pop.base$age==80],
                   b = pop.base$population[pop.base$sex==2&pop.base$age==84],
                   num.between = 4,
                   extrap.num = 1)
# Calculate m for men and women
grad.func <- function (X, num.groups, c) {
  (X - num.groups*c)/(num.groups*(num.groups-1)/2) # NB: This is to find the sum of 1,2,...,num.groups-1, since no gradient for first group.
}
m.m <- grad.func(X = pop.base$population[pop.base$sex==1&pop.base$age==85],
                 num.groups = 8, # Change this if the age groups change
                 c = c.m)
m.f <- grad.func(X = pop.base$population[pop.base$sex==2&pop.base$age==85],
                 num.groups = 10,
                 c = c.f)

new.entries <- NULL
# Generate new entries
for (i in 85:92) { # This should be the age groups that you want to extrapolate to
  Zm <- c(i,
          1,
          yr,
          m.m*(i-85)+c.m)
  new.entries <- rbind(new.entries,Zm)
}
for (i in 93:95){
  Zm <- c(i,
          1,
          yr,
          0)
  new.entries <- rbind(new.entries,Zm)
}
for (i in 85:94){
  Zf <- c(i,
          2,
          yr,
          m.f*(i-85)+c.f)
  new.entries <- rbind(new.entries,Zf)
}
for (i in 95:95){
  Zf <- c(i,
          2,
          yr,
          0)
  new.entries <- rbind(new.entries,Zf)
}

new.names <- c("age", "sex", "year", "population")
new.entries <- as.data.frame(new.entries)
colnames(new.entries) <- new.names

pop.base <- rbind(pop.base,new.entries)

# Plot the extended population database
par(mfrow=c(2,1))
plot(pop.base$age[pop.base$sex==1],pop.base$population[pop.base$sex==1],
     xlab = "Age", ylab = "Population, men")
plot(pop.base$age[pop.base$sex==2],pop.base$population[pop.base$sex==2],
     xlab = "Age", ylab = "Population, women")

# Finally remove the redundant point for 85+
pop.base <- pop.base[-c(86,172),]

#### Alternative code ####
# For the sake of estimating the median RR at each age point, it is not necessary to have an age distribution
# that matches the true age distribution - just necessary to have data points in each age category.
# To generate a dataset for this purpose, use code below.

#pop.base <- as.data.frame(cbind(rep(0:100,2),
#                                c(rep(1,101),rep(2,101)),
#                                rep(yr,202),
#                                rep(10,202)))
#
#new.names <- c("age", "sex", "year", "population")
#colnames(pop.base) <- new.names

######## Alternative code ends here
# Estimate this linear trend. (See Excel sheet, and load in new dataset) 
#if (work=="desktop") { # This sets the directory according to the initial settings
#  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\agent.fill")
#} else {
#  if (work=="laptop"){
#    setwd("F:\\microPRIME\\R modules\\agent.fill")
#  }  else { "Error" }
#}
#write.csv(pop.base, file = "population1998.csv") # THIS BIT NEEDS TO BE DONE IN EXCEL IF BASELINE YEAR CHANGED
#
# Step 3: Add in under 18s according to the fertility section of demo.R
#pop18     <- population[population$age==18&population$year>yr,2:5]
#pop18$age <- (yr+18) - pop18$year
## Merge this with the pop.base database
#pop18     <- pop18[order(pop18$age),]
#pop18$year<- yr
#write.csv(pop18, file = "population18_04OCT2019.csv") 
#
#
## Import extended baseline dataset
#pop.base <- read.csv("population1998extended.csv") # THIS BIT NEEDS TO BE DONE IN EXCEL IF BASELINE YEAR CHANGED

# Step 4: Turn into proportions of population
popwt.men        <- prop.table(pop.base$population[pop.base$sex==1&pop.base$age>=min.age]) 
popwt.women      <- prop.table(pop.base$population[pop.base$sex==2&pop.base$age>=min.age])
age.lookup.men   <- cbind(min.age:max.age, c(0, cumsum(popwt.men[1:(max.age-min.age)]))) # These generate a cumulative proportion of the population by age
age.lookup.women <- cbind(min.age:max.age, c(0, cumsum(popwt.women[1:(max.age-min.age)])))


## Generate populations with joint distributions of risk factors

# Start a loop around the m sets of results needed
# NB: Previous MI is NOT included, as this is added afterwards in accordance with logistic regression models
# calculated in risk.factors.R
for (i in 1:m){

### MEN
# Generate jointly distributed continuous variables
cop <- normalCopula(c(vector.sample$age.bmi.mr.m[i],  
                      vector.sample$age.sbp.mr.m[i],  
                      vector.sample$age.tcho.mr.m[i],
                      vector.sample$age.smok.mr.m[i],
                      vector.sample$age.diab.mr.m[i],
                      #vector.sample$age.mi.mr.m[i],
                      vector.sample$bmi.sbp.mr.m[i], 
                      vector.sample$bmi.tcho.mr.m[i],
                      vector.sample$bmi.smok.mr.m[i],
                      vector.sample$bmi.diab.mr.m[i],
                      #vector.sample$bmi.mi.mr.m[i],
                      vector.sample$sbp.tcho.mr.m[i],
                      vector.sample$sbp.smok.mr.m[i],
                      vector.sample$sbp.diab.mr.m[i],
                      #vector.sample$sbp.mi.mr.m[i],
                      vector.sample$tcho.smok.mr.m[i],
                      vector.sample$tcho.diab.mr.m[i],
                      #vector.sample$tcho.mi.mr.m[i],
                      vector.sample$smok.diab.mr.m[i]),
                      #vector.sample$smok.mi.mr.m[i],
                      #vector.sample$diab.mi.mr.m[i]),
                    dim = 6, dispstr="un")
x    <- rCopula(n, cop)

age.m  <- findInterval(x[,1],age.lookup.men[,2])
age.m  <- age.m + min.age - 1 # This is to translate the results so that the smallest return (which would
                              # be 1 untranslated) corresponds with the minimum age set earlier.

bmi.m  <- qlnorm(x[,2], vector.sample$bmi.pos.m[i],  vector.sample$bmi.shp.m[i]) # NB: update these transformations
sbp.m  <- qlnorm(x[,3], vector.sample$sbp.pos.m[i],  vector.sample$sbp.shp.m[i]) # if the distribuitions for the baseline
tcho.m <- qgamma(x[,4], vector.sample$tcho.pos.m[i], vector.sample$tcho.shp.m[i])# year in risk factors.R are different.

# For smoking and diabetes, the results are left as the [0,1] range that is correlated with the other variables. It is only
# the ranking that is important. In microsim, the trend by age and year will define the cutoff of people who smoke or have
# diabetes, but the rank will stay the same for each agent throughout.

#### OLD CODE FOR TRANSFORMING RANK TO SMOKING PROBABILITIES #######################################################
#smok.m <- x[,5]*vector.sample$smok.prev.m[i]/0.5                                                                  #
#diab.m <- x[,6]*vector.sample$diab.prev.m[i]/0.5                                                                  #
# NB: Above transformation only works because prevalence of smoking and diabetes is below 0.5 (therefore produces  #
# nothing outside of [0,1]). Note that this leaves the probability distribution as uniform, which is not ideal.    #
####################################################################################################################

smok.m <- x[,5]
diab.m <- x[,6]

#### OLD CODE FOR TRANSFORMING RANK TO MI PREVALENCE ###############################################################
# For previous MI, we just need a one or zero for whether or not there is a previous MI                            #
# mi.m <- c(rep(0,n)) # Holding                                                                                    #
#                                                                                                                  #     
#for (j in 1:n) {                                                # This translates the probability of MI           #
#  if (x[j,7]>(1-vector.sample$mi.prev.m[i])) {mi.m[j] <- 1}     # (effectively thrown out by the copula results)  #
#  else {mi.m[j] <- 0}                                           # into binary, MI or not.                         #
#}                                                                                                                 #
####################################################################################################################

#mi.m <- x[,7]

### WOMEN
# Generate jointly distributed continuous variables
cop <- normalCopula(c(vector.sample$age.bmi.mr.f[i],  
                      vector.sample$age.sbp.mr.f[i],  
                      vector.sample$age.tcho.mr.f[i],
                      vector.sample$age.smok.mr.f[i],
                      vector.sample$age.diab.mr.f[i],
                      #vector.sample$age.mi.mr.f[i],
                      vector.sample$bmi.sbp.mr.f[i], 
                      vector.sample$bmi.tcho.mr.f[i],
                      vector.sample$bmi.smok.mr.f[i],
                      vector.sample$bmi.diab.mr.f[i],
                      #vector.sample$bmi.mi.mr.f[i],
                      vector.sample$sbp.tcho.mr.f[i],
                      vector.sample$sbp.smok.mr.f[i],
                      vector.sample$sbp.diab.mr.f[i],
                      #vector.sample$sbp.mi.mr.f[i],
                      vector.sample$tcho.smok.mr.f[i],
                      vector.sample$tcho.diab.mr.f[i],
                      #vector.sample$tcho.mi.mr.f[i],
                      vector.sample$smok.diab.mr.f[i]),
                      #vector.sample$smok.mi.mr.f[i],
                      #vector.sample$diab.mi.mr.f[i]),
                    dim = 6, dispstr="un")
x    <- rCopula(n, cop)

age.f  <- findInterval(x[,1],age.lookup.women[,2])
age.f  <- age.f + min.age - 1 

bmi.f  <- qlnorm(x[,2], vector.sample$bmi.pos.f[i],  vector.sample$bmi.shp.f[i]) # NB: update these transformations
sbp.f  <- qlnorm(x[,3], vector.sample$sbp.pos.f[i],  vector.sample$sbp.shp.f[i]) # if the distribuitions for the baseline
tcho.f <- qgamma(x[,4], vector.sample$tcho.pos.f[i], vector.sample$tcho.shp.f[i])# year in risk factors.R are different.

smok.f <- x[,5]
diab.f <- x[,6]
#mi.f   <- x[,7]


# Export as csv
id.m               <- 1:n                # These few lines set up the IDs and sex assignation
male               <- c(rep(2,n))   # for the combined dataset to be used in microsim.
id.f               <- (n+1):(2*n)
female             <- c(rep(1,n)) # NB: This is where the categorisaton of sex comes in.
export.m           <- data.frame(cbind(id.m, 
                                       male, 
                                       age.m, 
                                       bmi.m, 
                                       sbp.m, 
                                       tcho.m, 
                                       smok.m, 
                                       diab.m))
                                       #mi.m))
export.f           <- data.frame(cbind(id.f, 
                                       female, 
                                       age.f, 
                                       bmi.f, 
                                       sbp.f, 
                                       tcho.f, 
                                       smok.f, 
                                       diab.f))
                                       #mi.f))
names              <- c("id", "sex2", "age", "bmi", "sbp", "tcho", "smok.p", "diab.p")# "mi")
colnames(export.m) <- names
colnames(export.f) <- names
export             <- merge(export.m, 
                            export.f, by = names, all = TRUE) # NB: these are not sorted on id
# Set directory for export
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Microsim inputs")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Microsim inputs")
  }  else { "Error" }
}

write.csv(export, paste0("input",i,"_best estimates_13MAR2020_18plus.csv")) # NB: Needs updating each time

# End loop
}

# model.run.end <- Sys.time() # This records the time when the loop was completed so I can calculate how long it took.

# model.run.end - model.run.start

# Add in the MI prevalence estimates by applying parameters from logistic regressions
check  <- read.csv("input1_best estimates_13MAR2020_18plus.csv") # Checking one set of outcomes

# Import the logistic regression results
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}

coefficients <- read.csv("miprev_coefficients_28NOV2019.csv")

expit <- function(x) {
  exp(x)/(1+exp(x))
}

check$mi2.p <- expit(coefficients[1,2]+
                       coefficients[2,2]*check$age+
                       coefficients[3,2]*check$age^2)
check$mi2   <- mapply(rbinom, 1,1,check$mi2.p)
check$mi1.p <- expit(coefficients[1,3]+
                       coefficients[2,3]*check$age+
                       coefficients[3,3]*check$age^2)
check$mi1   <- mapply(rbinom, 1,1,check$mi1.p)

sexchoose <- function(a,b,c) {
  if (a==1) { return(b) }
  else { return(c) }
}

check$mi <- mapply(sexchoose, check$sex2, check$mi1, check$mi2)

# Check that this looks appropriate
# Check prevalence of MI by sex and age group for comparison with external dataset
# Men
summary(check$mi[check$sex2==2&check$age>=55&check$age<=64])
summary(check$mi[check$sex2==2&check$age>=65&check$age<=74])
summary(check$mi[check$sex2==2&check$age>=75&check$age<=84])
summary(check$mi[check$sex2==2&check$age>=85])
# Women
summary(check$mi[check$sex2==1&check$age>=55&check$age<=64])
summary(check$mi[check$sex2==1&check$age>=65&check$age<=74])
summary(check$mi[check$sex2==1&check$age>=75&check$age<=84])
summary(check$mi[check$sex2==1&check$age>=85]) # These look pretty close

summary(check$age) # Works, min 41 max 95
test.m <- cbind(check$age[check$sex2==2],
                check$bmi[check$sex2==2],
                check$sbp[check$sex2==2],
                check$tcho[check$sex2==2],
                check$smok.p[check$sex2==2],
                check$diab.p[check$sex2==2],
                check$mi[check$sex2==2])
cor(test.m) # Reproduces similar marginal correlations

test.f <- cbind(check$age[check$sex==1],
                check$bmi[check$sex==1],
                check$sbp[check$sex==1],
                check$tcho[check$sex==1],
                check$smok.p[check$sex==1],
                check$diab.p[check$sex==1],
                check$mi[check$sex==1])
cor(test.f) # Reproduces similar marginal correlations

summary(check$smok.p[check$sex2==2])
summary(check$smok.p[check$sex2==1])
summary(check$diab.p[check$sex2==2])
summary(check$diab.p[check$sex2==1])
summary(check$mi[check$sex2==2])
summary(check$mi[check$sex2==1])

# Overwrite the final dataset
export2 <- as.data.frame(cbind(check$id,
                               check$sex2,
                               check$age,
                               check$bmi,
                               check$sbp,
                               check$tcho,
                               check$smok.p,
                               check$diab.p,
                               check$mi))
names <- c("id","sex2","age","bmi","sbp","tcho","smok.p","diab.p","mi")
colnames(export2) <- names

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Microsim inputs")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Microsim inputs")
  }  else { "Error" }
}

write.csv(export2, "input1_best estimates_13MAR2020_18plus.csv") # NB: Needs updating each time

# Reduce to over 35s:
export3 <- export2[export2$age>=35,]

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Microsim inputs")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Microsim inputs")
  }  else { "Error" }
}

write.csv(export3, "input1_best estimates_13MAR2020_35plus.csv") # NB: Needs updating each time



### NOTES

# MIGRATION: This will be dealt with in the microsim model, by generating probabilities of new
# people, drawn from a separate sample of migrants, when necessary in each year.

# This R script has been updated from agent fill.R so that smoking and diabetes are just
# probabilities instead of actual binary estimates.