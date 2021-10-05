### risk factors ###

# Last updated: 25th February 2020

# This R script aims to generate parameters to be transferred into Parameter Space. This includes: 
# 1. estimates of mean, sd and 99% confidence intervals for distributions of characteristics
# 2. estimates of prevalence and 99% confidence intervals of binary risk factors
# 3. joint distributions of risk factors for agent fill.R
# 4. descriptions of trends in risk factors by age and year, stratified by sex

# Install packages
# install.packages("fitdistrplus")
# install.packages("weights")
# install.packages("base")
# install.packages("haven")
# install.packages("plyr")
# install.packages("lmtest")
library("fitdistrplus")
library("weights")
library("Hmisc")
library("base")
library("haven")
library("plyr")
library("lmtest")

# Set work environment
work <- "desktop" # Choice between 'desktop' and 'laptop'. Ensures script uses correct directories.

##################################################################################################
#                                                                                                #
# 1. Estimates of distributions of continuous risk factors in baseline year (1992)               #
#                                                                                                #
##################################################################################################

# Import HSE master dataset
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Health Survey for England")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Health Survey for England")
  }  else { "Error" }
}
hse.master   <- read_dta('hse_master.dta')
hse.master$x <- 1                 # This is needed for next line
n            <- sum(hse.master$x) # This sets n as the number of lines in hse.master

# Generate complete year variable
yearconvert <- function(x){ 
  if (x<90){ 2000 + x }           
  else { 1900 + x }
}

hse.master$year2 <- lapply(hse.master$year, yearconvert)
hse.master$year2 <- as.numeric(hse.master$year2)

# Check data in hse.master
# Data points by year
table(hse.master$x, hse.master$year2) # Some years have many more data points. Boosts in those years?
# Distribution of age in whole database
summary(hse.master$age) # All aged 18 and over.
hist(hse.master$age)
# Smoking data by year
summary(hse.master$cig1)
table(hse.master$x[is.na(hse.master$cig1)==FALSE], hse.master$year2[is.na(hse.master$cig1)==FALSE])
# Diabetes data by year
summary(hse.master$diabetes) 
table(hse.master$x[is.na(hse.master$diabetes)==FALSE], hse.master$year2[is.na(hse.master$diabetes)==FALSE])
# BMI data by year
summary(hse.master$bmi)
table(hse.master$x[is.na(hse.master$bmi)==FALSE], hse.master$year2[is.na(hse.master$bmi)==FALSE])
# SBP data by year
summary(hse.master$sbp)
table(hse.master$x[is.na(hse.master$sbp)==FALSE], hse.master$year2[is.na(hse.master$sbp)==FALSE])
# Cholesterol data by year
summary(hse.master$chol) 
table(hse.master$x[is.na(hse.master$chol)==FALSE], hse.master$year2[is.na(hse.master$chol)==FALSE])

# Set year for data extraction
yr <- 2006 # For extraction of the initial distributions of the continuous risk factors and the prevalence of the 
           # binary risk factors, this will be based on data from three years centred at yr. 


# THIS NEXT SECTION IS DESIGNED AS A CHECK THAT WEIGHTING THE DATA BY SURVEY WEIGHTS AND BY 
# POPULATION WEIGHTS (i.e. age-standardising to the England population) DOES NOT AFFECT
# DISTRIBUTIONS OF CONTINUOUS VARIABLES, AS THE FITDISTPLUS DOES NOT WORK WITH WEIGHTING.
##############################################################################################

# Import the ONS population estimates
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Demographic trends")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Demographic trends")
  }  else { "Error" }
}

fert.mort <- read.csv("pop.mort04AUG2017.csv") # To update the input file, need to run 'demo_preparing dataset.R'

# Create dataset of over 18s for identified yr
eng.pop <- fert.mort[fert.mort$year==yr & fert.mort$age>=18,]
par(mfrow=c(2,1)) # Plots of English population by age
barplot(eng.pop$population[eng.pop$sex==1], xlab = "age", ylab = "pop", main = "Men")
barplot(eng.pop$population[eng.pop$sex==2], xlab = "age", ylab = "pop", main = "Women")

# Convert to proportions of over 18s
popwt.men   <- prop.table(eng.pop$population[eng.pop$sex==1])
popwt.women <- prop.table(eng.pop$population[eng.pop$sex==2])

# Get count of people in HSfE
eightyfiveplus <- function(x)   { # This will convert over 85 ages in HSfE to be comparable with
  if (is.na(x)==TRUE) { NA }
  else if (x<85) { x }              # English  population  estimates
  else { 85 }
}

hse.master$age2 <- lapply(hse.master$age, eightyfiveplus)

hse.master$age2 <- 1 # Holding variable
for (i in 1:n){
  hse.master$age2[i] <- eightyfiveplus(hse.master$age[i], y) 
}

hsfe.men          <- prop.table(table(hse.master$age2[hse.master$sex==1&hse.master$year2==yr]))
hsfe.women        <- prop.table(table(hse.master$age2[hse.master$sex==2&hse.master$year2==yr]))
wt.men            <- popwt.men   / hsfe.men
wt.women          <- popwt.women / hsfe.women # Weight variables by age2 in hse.master
weights           <- data.frame(cbind(18:85, as.numeric(wt.men), as.numeric(wt.women)))
names             <- c("age2", "wt.men", "wt.women")
colnames(weights) <- names

# Merge dataset
hseyr <- hse.master[hse.master$year2==yr,]
hse.merge <- merge(hseyr, weights, by = 'age2')

# Checking distribution by weighting: Unweighted versus weighted histogram of BMI, men 
par(mfcol = c(2,1)) # MEN
hist(hse.merge$bmi[hse.merge$sex==1],     breaks = 25, main = "Unweighted")
wtd.hist(hse.merge$bmi[hse.merge$sex==1], breaks = 25, main = "Population weights", 
         weight = hse.merge$wt.men[hse.merge$sex==1])

par(mfcol = c(2,1)) # WOMEN
hist(hse.merge$bmi[hse.merge$sex==2],     breaks = 25, main = "Unweighted")
wtd.hist(hse.merge$bmi[hse.merge$sex==2], breaks = 25, main = "Population weights", 
         weight = hse.merge$wt.women[hse.merge$sex==2])

# DISTRIBUTIONS SHOW THAT THERE IS VERY LITTLE DIFFERENCE IN THE DISTRIBUTION OF BMI
# WHETHER UNWEIGHTED, WEIGHTED BY SURVEY WEIGHTS, OR WEIGHTED BY THE ENGLISH  POPULATION IN 2015
# THIS MEANS IT IS SAFE TO PROCEED WITH UNWEIGHTED SURVEY DATA FOR CALCULATING DISTRIBUTIONS
# OF CONTINUOUS RISK FACTORS
# 
#############################################################################################


# Return to calculations of shape and position parameters

par(mfcol = c(1,1)) # Reset the graph outputs to 1 per page.

# BMI: MEN
# Check distribution
bmi.men <- as.numeric(na.omit(hse.master$bmi[hse.master$sex==1                 # Sex is male    AND
                                             & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                | hse.master$year2==yr         # Year is yr     OR
                                                | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(bmi.men) # Plot suggests either lognormal, gamma or weibull
lnbmi.men      <- fitdist(bmi.men, "lnorm")
gammabmi.men   <- fitdist(bmi.men, "gamma")
weibullbmi.men <- fitdist(bmi.men, "weibull")
gofstat(list(lnbmi.men, gammabmi.men, weibullbmi.men), fitnames = c("lognormal", "gamma", "weibull"))
# Lognormal has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
bmi.pos.m <- c(as.numeric(lnbmi.men$estimate[1]-2.58*lnbmi.men$sd[1]),
               as.numeric(lnbmi.men$estimate[1]+2.58*lnbmi.men$sd[1]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")
bmi.shp.m <- c(as.numeric(lnbmi.men$estimate[2]-2.58*lnbmi.men$sd[2]),
               as.numeric(lnbmi.men$estimate[2]+2.58*lnbmi.men$sd[2]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")

# BMI: WOMEN
# Check distribution
bmi.women <- as.numeric(na.omit(hse.master$bmi[hse.master$sex==2                 # Sex is female  AND
                                               & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                  | hse.master$year2==yr         # Year is yr     OR
                                                  | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(bmi.women) # Plot suggests either beta, gamma, lognormal or weibull NB: Beta should be between [0,1]
lnbmi.women      <- fitdist(bmi.women, "lnorm")
gammabmi.women   <- fitdist(bmi.women, "gamma")
weibullbmi.women <- fitdist(bmi.women, "weibull")
gofstat(list(lnbmi.women, gammabmi.women, weibullbmi.women), 
        fitnames = c("lognormal", "gamma", "weibull"))
# Lognormal has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
bmi.pos.f <- c(as.numeric(lnbmi.women$estimate[1]-2.58*lnbmi.women$sd[1]),
               as.numeric(lnbmi.women$estimate[1]+2.58*lnbmi.women$sd[1]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")
bmi.shp.f <- c(as.numeric(lnbmi.women$estimate[2]-2.58*lnbmi.women$sd[2]),
               as.numeric(lnbmi.women$estimate[2]+2.58*lnbmi.women$sd[2]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")

# SBP: MEN
# Check distribution
sbp.men <- as.numeric(na.omit(hse.master$sbp[hse.master$sex==1                 # Sex is male    AND
                                             & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                | hse.master$year2==yr         # Year is yr     OR
                                                | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(sbp.men) # Plot suggests either lognormal, gamma or weibull. Try normal as well

nsbp.men       <- fitdist(sbp.men, "norm")
lnsbp.men      <- fitdist(sbp.men, "lnorm")
gammasbp.men   <- fitdist(sbp.men, "gamma")
weibullsbp.men <- fitdist(sbp.men, "weibull")
gofstat(list(nsbp.men, lnsbp.men, gammasbp.men, weibullsbp.men), 
        fitnames = c("normal", "lognormal", "gamma", "weibull"))
# Lognormal has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
sbp.pos.m <- c(as.numeric(lnsbp.men$estimate[1]-2.58*lnsbp.men$sd[1]),
               as.numeric(lnsbp.men$estimate[1]+2.58*lnsbp.men$sd[1]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")
sbp.shp.m <- c(as.numeric(lnsbp.men$estimate[2]-2.58*lnsbp.men$sd[2]),
               as.numeric(lnsbp.men$estimate[2]+2.58*lnsbp.men$sd[2]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")

# SBP: WOMEN
# Check distribution
sbp.women <- as.numeric(na.omit(hse.master$sbp[hse.master$sex==2                 # Sex is female  AND
                                               & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                  | hse.master$year2==yr         # Year is yr     OR
                                                  | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(sbp.women) # Plot suggests either lognormal, gamma or weibull. Try normal as well

nsbp.women       <- fitdist(sbp.women, "norm")
lnsbp.women      <- fitdist(sbp.women, "lnorm")
gammasbp.women   <- fitdist(sbp.women, "gamma")
weibullsbp.women <- fitdist(sbp.women, "weibull")
gofstat(list(nsbp.women, lnsbp.women, gammasbp.women, weibullsbp.women), 
        fitnames = c("normal", "lognormal", "gamma", "weibull"))
# Lognormal has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
sbp.pos.f <- c(as.numeric(lnsbp.women$estimate[1]-2.58*lnsbp.women$sd[1]),
               as.numeric(lnsbp.women$estimate[1]+2.58*lnsbp.women$sd[1]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")
sbp.shp.f <- c(as.numeric(lnsbp.women$estimate[2]-2.58*lnsbp.women$sd[2]),
               as.numeric(lnsbp.women$estimate[2]+2.58*lnsbp.women$sd[2]),
               NA, NA, NA, NA, NA, NA, NA, NA, "lnorm")

# CHOLESTEROL: MEN
# Check distribution
tcho.men <- as.numeric(na.omit(hse.master$chol[hse.master$sex==1                 # Sex is male    AND
                                               & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                  | hse.master$year2==yr         # Year is yr     OR
                                                  | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(tcho.men) # Plot suggests normal. Try lognormal, gamma or weibull as well

ntcho.men       <- fitdist(tcho.men, "norm")
lntcho.men      <- fitdist(tcho.men, "lnorm")
gammatcho.men   <- fitdist(tcho.men, "gamma")
weibulltcho.men <- fitdist(tcho.men, "weibull")
gofstat(list(ntcho.men, lntcho.men, gammatcho.men, weibulltcho.men), 
        fitnames = c("normal", "lognormal", "gamma", "weibull"))
# Gamma has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
tcho.pos.m <- c(as.numeric(gammatcho.men$estimate[1]-2.58*gammatcho.men$sd[1]),
               as.numeric(gammatcho.men$estimate[1]+2.58*gammatcho.men$sd[1]),
               NA, NA, NA, NA, NA, NA, NA, NA, "gamma")
tcho.shp.m <- c(as.numeric(gammatcho.men$estimate[2]-2.58*gammatcho.men$sd[2]),
               as.numeric(gammatcho.men$estimate[2]+2.58*gammatcho.men$sd[2]),
               NA, NA, NA, NA, NA, NA, NA, NA, "gamma")

# CHOLESTEROL: WOMEN
# Check distribution
tcho.women <- as.numeric(na.omit(hse.master$chol[hse.master$sex==2                 # Sex is female  AND
                                                 & (hse.master$year2==(yr-1)       # Year is (yr-1) OR
                                                    | hse.master$year2==yr         # Year is yr     OR
                                                    | hse.master$year2==(yr+1))])) # Year is (yr+1)
descdist(tcho.women) # Plot suggests gamma. Try lognormal, gamma, normal or weibull as well

ntcho.women       <- fitdist(tcho.women, "norm")
lntcho.women      <- fitdist(tcho.women, "lnorm")
gammatcho.women   <- fitdist(tcho.women, "gamma")
weibulltcho.women <- fitdist(tcho.women, "weibull")
gofstat(list(ntcho.women, lntcho.women, gammatcho.women, weibulltcho.women), 
        fitnames = c("normal", "lognormal", "gamma", "weibull"))
# Gamma has lowest Kologorov-Smirnov statistic.

## Prepare data for exporting
# NB: 99% confidence interval for mean and and sd are based on mean +/- 2.58 * Standard error
tcho.pos.f <- c(as.numeric(gammatcho.women$estimate[1]-2.58*gammatcho.women$sd[1]),
                as.numeric(gammatcho.women$estimate[1]+2.58*gammatcho.women$sd[1]),
                NA, NA, NA, NA, NA, NA, NA, NA, "gamma")
tcho.shp.f <- c(as.numeric(gammatcho.women$estimate[2]-2.58*gammatcho.women$sd[2]),
                as.numeric(gammatcho.women$estimate[2]+2.58*gammatcho.women$sd[2]),
                NA, NA, NA, NA, NA, NA, NA, NA, "gamma")

##################################################################################################
#                                                                                                #
# 2. Prevalence of binary risk factors in baseline year                                          #
#                                                                                                #
##################################################################################################

# Generate new smoking variable {smok: 0 = not current smoker; 1= current smoker}
smokconvert <- function(x){
  if (is.na(x)==TRUE) { NA }
  else if (x==3) { 1 }
    else { 0 }
}
  

hse.master$smok <- as.numeric(lapply(hse.master$cig1, smokconvert))

# Generate new previous MI variable {mi: 0 = no previous MI; 1 = previous MI}
miconvert <- function(x){
  if (is.na(x)==TRUE) { NA }
  else if (x==1|x==3) { 1 }
  else { 0 }
}

hse.master$mi <- as.numeric(lapply(hse.master$doc_diagnosedMI_ever, miconvert))

# Generate proportions of smokers by five year age group and sex in three year circle round yr
hse3yr    <- data.frame(hse.master[hse.master$year2==(yr-1)|hse.master$year2==yr|hse.master$year2==(yr+1),])

# NB: Section below produces age-standardised prevalence of smoking, diabetes and MI prevalence.
# Removed for model 11, as for copula you just need the unadjusted prevalence in population (this is
# because the copula includes joint distributions with age)
# Age-adjusted prevalence of smoking
# NB: age weighting is by 5 years until age 90, and then all pop.
# Age weights sum to 1 (each is just number of people divided by total population in yr)
#
#ageweights   <- c("18-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69",
#                  "70-74","75-79","80-84","85+")
#menweights   <- c(sum(eng.pop[1:7,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[8:12,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[13:17,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[18:22,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[23:27,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[28:32,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[33:37,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[38:42,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[43:47,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[48:52,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[53:57,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[58:62,5])/sum(eng.pop[1:68,5]),
#                  sum(eng.pop[63:67,5])/sum(eng.pop[1:68,5]),
#                  eng.pop[68,5]/sum(eng.pop[1:68,5]))
#womenweights <- c(sum(eng.pop[69:75,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[76:80,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[81:85,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[86:90,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[91:95,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[96:100,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[101:105,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[106:110,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[111:115,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[116:120,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[121:125,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[126:130,5])/sum(eng.pop[69:136,5]),
#                  sum(eng.pop[131:135,5])/sum(eng.pop[69:136,5]),
#                  eng.pop[136,5]/sum(eng.pop[69:136,5]))
#five.yr.wt <- cbind(ageweights, menweights, womenweights)
#
#myprop <- function(x, y){ # x is the variable needed for age-sex-specific proportions, y = 1 or 2 for the sex
#  c(sum(na.omit(x[hse3yr$sex==y&hse3yr$age<=24]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age<=24])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=25&hse3yr$age<=29]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=25&hse3yr$age<=29])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=30&hse3yr$age<=34]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=30&hse3yr$age<=34])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=35&hse3yr$age<=39]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=35&hse3yr$age<=39])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=40&hse3yr$age<=44]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=40&hse3yr$age<=44])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=45&hse3yr$age<=49]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=45&hse3yr$age<=49])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=50&hse3yr$age<=54]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=50&hse3yr$age<=54])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=55&hse3yr$age<=59]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=55&hse3yr$age<=59])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=60&hse3yr$age<=64]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=60&hse3yr$age<=64])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=65&hse3yr$age<=69]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=65&hse3yr$age<=69])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=70&hse3yr$age<=74]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=70&hse3yr$age<=74])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=75&hse3yr$age<=79]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=75&hse3yr$age<=79])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=80&hse3yr$age<=84]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=80&hse3yr$age<=84])),
#    sum(na.omit(x[hse3yr$sex==y&hse3yr$age>=85]))/sum(!is.na(x[hse3yr$sex==y&hse3yr$age>=85])))
#}
#
#mensmok   <- myprop(hse3yr$smok, 1)
#womensmok <- myprop(hse3yr$smok, 2)
#smok.prev <- cbind(ageweights, mensmok, womensmok)

## Prepare data for exporting
#smok.prev.m    <- sum(as.numeric(smok.prev[,2])*as.numeric(five.yr.wt[,2]))
#smok.prev.m.se <- sqrt((smok.prev.m*(1-smok.prev.m)*sum(as.numeric(five.yr.wt[,2])^2))/sum(!is.na(hse3yr$smok[hse3yr$sex==1])))
#smok.prev.m    <- c(smok.prev.m-2.58*smok.prev.m.se, smok.prev.m+2.58*smok.prev.m.se,
#                    NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
#smok.prev.f    <- sum(as.numeric(smok.prev[,3])*as.numeric(five.yr.wt[,3]))
#smok.prev.f.se <- sqrt((smok.prev.f*(1-smok.prev.f)*sum(as.numeric(five.yr.wt[,3])^2))/sum(!is.na(hse3yr$smok[hse3yr$sex==2])))
#smok.prev.f    <- c(smok.prev.f-2.58*smok.prev.f.se, smok.prev.f+2.58*smok.prev.f.se,
#                    NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
#
## Diabetes prevalence
#mendiab   <- myprop(hse3yr$diabetes, 1)
#womendiab <- myprop(hse3yr$diabetes, 2)
#diab.prev <- cbind(ageweights, mendiab, womendiab) # Note that with 1992 data there is not much of a relationship between
#                                                   # age and diabetes.
#
## Prepare data for exporting
#diab.prev.m    <- sum(as.numeric(diab.prev[,2])*as.numeric(five.yr.wt[,2]))
#diab.prev.m.se <- sqrt((diab.prev.m*(1-diab.prev.m)*sum(as.numeric(five.yr.wt[,2])^2))/sum(!is.na(hse3yr$diabetes[hse3yr$sex==1])))
#diab.prev.m    <- c(diab.prev.m-2.58*diab.prev.m.se, diab.prev.m+2.58*diab.prev.m.se,
#                    NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
#diab.prev.f    <- sum(as.numeric(diab.prev[,3])*as.numeric(five.yr.wt[,3]))
#diab.prev.f.se <- sqrt((diab.prev.f*(1-diab.prev.f)*sum(as.numeric(five.yr.wt[,3])^2))/sum(!is.na(hse3yr$diabetes[hse3yr$sex==2])))
#diab.prev.f    <- c(diab.prev.f-2.58*diab.prev.f.se, diab.prev.f+2.58*diab.prev.f.se,
#                    NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
#
## Prevalence of previous heart attack
#menprevmi   <- myprop(hse3yr$mi, 1)
#womenprevmi <- myprop(hse3yr$mi, 2)
#mi.prev <- cbind(ageweights, menprevmi, womenprevmi)
#
## Prepare data for exporting
#mi.prev.m    <- sum(as.numeric(mi.prev[,2])*as.numeric(five.yr.wt[,2]))
#mi.prev.m.se <- sqrt((mi.prev.m*(1-mi.prev.m)*sum(as.numeric(five.yr.wt[,2])^2))/sum(!is.na(hse3yr$mi[hse3yr$sex==1])))
#mi.prev.m    <- c(mi.prev.m-2.58*mi.prev.m.se, mi.prev.m+2.58*mi.prev.m.se,
#                  NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
#mi.prev.f    <- sum(as.numeric(mi.prev[,3])*as.numeric(five.yr.wt[,3]))
#mi.prev.f.se <- sqrt((mi.prev.f*(1-mi.prev.f)*sum(as.numeric(five.yr.wt[,3])^2))/sum(!is.na(hse3yr$mi[hse3yr$sex==2])))
#mi.prev.f    <- c(mi.prev.f-2.58*mi.prev.f.se, mi.prev.f+2.58*mi.prev.f.se,
#                  NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")

# Smoking prevalence data
smok.prev.m   <- sum(hse3yr$x[hse3yr$sex==1&hse3yr$smok==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE)
smok.prev.m.se<- sqrt(smok.prev.m*(1-smok.prev.m)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE))
smok.prev.m   <- c(smok.prev.m-2.58*smok.prev.m.se, smok.prev.m+2.58*smok.prev.m.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
smok.prev.f   <- sum(hse3yr$x[hse3yr$sex==2&hse3yr$smok==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE)
smok.prev.f.se<- sqrt(smok.prev.f*(1-smok.prev.f)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE))
smok.prev.f   <- c(smok.prev.f-2.58*smok.prev.f.se, smok.prev.f+2.58*smok.prev.f.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")

# Diabetes prevalence data
diab.prev.m   <- sum(hse3yr$x[hse3yr$sex==1&hse3yr$diab==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE)
diab.prev.m.se<- sqrt(diab.prev.m*(1-diab.prev.m)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE))
diab.prev.m   <- c(diab.prev.m-2.58*diab.prev.m.se, diab.prev.m+2.58*diab.prev.m.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
diab.prev.f   <- sum(hse3yr$x[hse3yr$sex==2&hse3yr$diab==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE)
diab.prev.f.se<- sqrt(diab.prev.f*(1-diab.prev.f)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE))
diab.prev.f   <- c(diab.prev.f-2.58*diab.prev.f.se, diab.prev.f+2.58*diab.prev.f.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")

# MI prevalence data
mi.prev.m     <- sum(hse3yr$x[hse3yr$sex==1&hse3yr$mi==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE)
mi.prev.m.se  <- sqrt(mi.prev.m*(1-mi.prev.m)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE))
mi.prev.m     <- c(mi.prev.m-2.58*mi.prev.m.se, mi.prev.m+2.58*mi.prev.m.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")
mi.prev.f     <- sum(hse3yr$x[hse3yr$sex==2&hse3yr$mi==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE)
mi.prev.f.se  <- sqrt(mi.prev.f*(1-mi.prev.f)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE))
mi.prev.f     <- c(mi.prev.f-2.58*mi.prev.f.se, mi.prev.f+2.58*mi.prev.f.se,
                   NA, NA, NA, NA, NA, NA, NA, NA, "prevalence")

# Check age-specific prevalence of previous MI
# Men
summary(hse3yr$mi[hse3yr$sex==1&hse3yr$age>=55&hse3yr$age<=64])
summary(hse3yr$mi[hse3yr$sex==1&hse3yr$age>=65&hse3yr$age<=74])
summary(hse3yr$mi[hse3yr$sex==1&hse3yr$age>=75&hse3yr$age<=84])
summary(hse3yr$mi[hse3yr$sex==1&hse3yr$age>=85])
# Women
summary(hse3yr$mi[hse3yr$sex==2&hse3yr$age>=55&hse3yr$age<=64])
summary(hse3yr$mi[hse3yr$sex==2&hse3yr$age>=65&hse3yr$age<=74])
summary(hse3yr$mi[hse3yr$sex==2&hse3yr$age>=75&hse3yr$age<=84])
summary(hse3yr$mi[hse3yr$sex==2&hse3yr$age>=85])



##################################################################################################
#                                                                                                #
# 3. Joint distributions of risk factors                                                         #
#                                                                                                #
##################################################################################################

## Marginal correlations of continuous risk factors

# NB: See 'riSk factors_checking prevalence rates.R' and 'agent fill_checking prevalence rates.R' for
# in depth tests of different methods for generating parameters to ensure appropriate joint distributions
# and prevalence rates. The preferred method is as follows:
# - Age included in the correlations
# - Actual values used for continuous variables (i.e. not rank correlation)
# - Binary variables included in the correlation matrices (i.e. rather than logistic regression models
#   which faiiled to reproduce prevalence rates)
#
# The correlations are drawn from the three years around the baseline year. But the trend analyses will 
# not allow for the joint correlations to change over time. But different trends in the risk factors will
# result in changes in joint distribution over time anyway.

# Age:BMI marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$bmi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  # This stores the 99% confidence intervals 
                                                     # in a vector
age.bmi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")       # This prepares the results for export
age.bmi.mr.m.estimate <- as.numeric(A$estimate)      # This stores the 'best estimate'

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$bmi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.bmi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")       # This prepares the results for export
age.bmi.mr.f.estimate <- as.numeric(A$estimate)


# Age:SBP marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$sbp[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.sbp.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.sbp.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$sbp[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.sbp.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.sbp.mr.f.estimate <- as.numeric(A$estimate)

# Age:Total cholesterol marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$chol[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.tcho.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.tcho.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$chol[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.tcho.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.tcho.mr.f.estimate <- as.numeric(A$estimate)

# Age:Smoking marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$smok[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.smok.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.smok.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$smok[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.smok.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.smok.mr.f.estimate <- as.numeric(A$estimate)

# Age:Diabetes marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$diabetes[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.diab.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.diab.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$diabetes[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.diab.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
age.diab.mr.f.estimate <- as.numeric(A$estimate)

# Age:Previous MI marginal correlation
A <- cor.test(hse3yr$age[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
age.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$age[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
age.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
age.mi.mr.f.estimate <- as.numeric(A$estimate)


# BMI: SBP
A <- cor.test(hse3yr$bmi[hse3yr$sex==1],hse3yr$sbp[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99) # This does the correlation
B <- matrix(unlist(A[9]), ncol = 1)                  # This stores the 99% confidence intervals 
                                                     # in a vector
bmi.sbp.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")       # This prepares the results for export
bmi.sbp.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$bmi[hse3yr$sex==2],hse3yr$sbp[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)
bmi.sbp.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, NA, NA, NA, NA, "pearson r")
bmi.sbp.mr.f.estimate <- as.numeric(A$estimate)

# BMI: Total cholesterol
A <- cor.test(hse3yr$bmi[hse3yr$sex==1],hse3yr$chol[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)
bmi.tcho.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, NA, NA, NA, NA, "pearson r")
bmi.tcho.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$bmi[hse3yr$sex==2],hse3yr$chol[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)
bmi.tcho.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, NA, NA, NA, NA, "pearson r")
bmi.tcho.mr.f.estimate <- as.numeric(A$estimate)

# BMI: Smoking marginal correlation
A <- cor.test(hse3yr$bmi[hse3yr$sex==1],hse3yr$smok[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.smok.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
bmi.smok.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$bmi[hse3yr$sex==2],hse3yr$smok[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.smok.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                  NA, NA, NA, NA, "pearson r")
bmi.smok.mr.f.estimate <- as.numeric(A$estimate)

# BMI: Diabetes marginal correlation
A <- cor.test(hse3yr$bmi[hse3yr$sex==1],hse3yr$diabetes[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.diab.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
bmi.diab.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$bmi[hse3yr$sex==2],hse3yr$diabetes[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.diab.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
bmi.diab.mr.f.estimate <- as.numeric(A$estimate)

# BMI: Previous MI marginal correlation
A <- cor.test(hse3yr$bmi[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
bmi.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$bmi[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
bmi.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
bmi.mi.mr.f.estimate <- as.numeric(A$estimate)

# SBP: Total cholesterol
A <- cor.test(hse3yr$sbp[hse3yr$sex==1],hse3yr$chol[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)
sbp.tcho.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, NA, NA, NA, NA, "pearson r")
sbp.tcho.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$sbp[hse3yr$sex==2],hse3yr$chol[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)
sbp.tcho.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, NA, NA, NA, NA, "pearson r")
sbp.tcho.mr.f.estimate <- as.numeric(A$estimate)

# SBP: Smoking marginal correlation
A <- cor.test(hse3yr$sbp[hse3yr$sex==1],hse3yr$smok[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.smok.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
sbp.smok.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$sbp[hse3yr$sex==2],hse3yr$smok[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.smok.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r") 
sbp.smok.mr.f.estimate <- as.numeric(A$estimate)

# SBP: Diabetes marginal correlation
A <- cor.test(hse3yr$sbp[hse3yr$sex==1],hse3yr$diabetes[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.diab.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
sbp.diab.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$sbp[hse3yr$sex==2],hse3yr$diabetes[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.diab.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")     
sbp.diab.mr.f.estimate <- as.numeric(A$estimate)

# SBP: Previous MI marginal correlation
A <- cor.test(hse3yr$sbp[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
sbp.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$sbp[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
sbp.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")     
sbp.mi.mr.f.estimate <- as.numeric(A$estimate)

# Total cholesterol: Smoking marginal correlation
A <- cor.test(hse3yr$chol[hse3yr$sex==1],hse3yr$smok[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.smok.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
tcho.smok.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$chol[hse3yr$sex==2],hse3yr$smok[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.smok.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
tcho.smok.mr.f.estimate <- as.numeric(A$estimate)

# Total cholesterol: Diabetes marginal correlation
A <- cor.test(hse3yr$chol[hse3yr$sex==1],hse3yr$diabetes[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.diab.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
tcho.diab.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$chol[hse3yr$sex==2],hse3yr$diabetes[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.diab.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")  
tcho.diab.mr.f.estimate <- as.numeric(A$estimate)

# Total cholesterol: Previous MI marginal correlation
A <- cor.test(hse3yr$chol[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")
tcho.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$chol[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
tcho.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")  
tcho.mi.mr.f.estimate <- as.numeric(A$estimate)

# Smoking: Diabetes marginal correlation
A <- cor.test(hse3yr$smok[hse3yr$sex==1],hse3yr$diabetes[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
smok.diab.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")
smok.diab.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$smok[hse3yr$sex==2],hse3yr$diabetes[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
smok.diab.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                   NA, NA, NA, NA, "pearson r")  
smok.diab.mr.f.estimate <- as.numeric(A$estimate)

# Smoking: Previous MI marginal correlation
A <- cor.test(hse3yr$smok[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
smok.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")
smok.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$smok[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
smok.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")  
smok.mi.mr.f.estimate <- as.numeric(A$estimate)

# Diabetes: Previous MI marginal correlation
A <- cor.test(hse3yr$diabetes[hse3yr$sex==1],hse3yr$mi[hse3yr$sex==1], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
diab.mi.mr.m <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")
diab.mi.mr.m.estimate <- as.numeric(A$estimate)

A <- cor.test(hse3yr$diabetes[hse3yr$sex==2],hse3yr$mi[hse3yr$sex==2], 
              method = "pearson", conf.level = 0.99)
B <- matrix(unlist(A[9]), ncol = 1)                  
diab.mi.mr.f <- c(B[1,1], B[2,1], NA, NA, NA, NA, 
                    NA, NA, NA, NA, "pearson r")  
diab.mi.mr.f.estimate <- as.numeric(A$estimate)


# Combine into single dataset
parameter.space <- rbind(bmi.pos.m,
                         bmi.shp.m,
                         bmi.pos.f,
                         bmi.shp.f,
                         sbp.pos.m,
                         sbp.shp.m,
                         sbp.pos.f,
                         sbp.shp.f,
                         tcho.pos.m,
                         tcho.shp.m,
                         tcho.pos.f,
                         tcho.shp.f,
                         smok.prev.m,
                         smok.prev.f,
                         diab.prev.m,
                         diab.prev.f,
                         mi.prev.m,
                         mi.prev.f,
                         age.bmi.mr.m,
                         age.sbp.mr.m,
                         age.tcho.mr.m,
                         age.smok.mr.m,
                         age.diab.mr.m,
                         age.mi.mr.m,
                         bmi.sbp.mr.m,
                         bmi.tcho.mr.m,
                         bmi.smok.mr.m,
                         bmi.diab.mr.m,
                         bmi.mi.mr.m,
                         sbp.tcho.mr.m,
                         sbp.smok.mr.m,
                         sbp.diab.mr.m,
                         sbp.mi.mr.m,
                         tcho.smok.mr.m,
                         tcho.diab.mr.m,
                         tcho.mi.mr.m,
                         smok.diab.mr.m,
                         smok.mi.mr.m,
                         diab.mi.mr.m,
                         age.bmi.mr.f,
                         age.sbp.mr.f,
                         age.tcho.mr.f,
                         age.smok.mr.f,
                         age.diab.mr.f,
                         age.mi.mr.f,
                         bmi.sbp.mr.f,
                         bmi.tcho.mr.f,
                         bmi.smok.mr.f,
                         bmi.diab.mr.f,
                         bmi.mi.mr.f,
                         sbp.tcho.mr.f,
                         sbp.smok.mr.f,
                         sbp.diab.mr.f,
                         sbp.mi.mr.f,
                         tcho.smok.mr.f,
                         tcho.diab.mr.f,
                         tcho.mi.mr.f,
                         smok.diab.mr.f,
                         smok.mi.mr.f,
                         diab.mi.mr.f)
varnames <- c("range1lo",
              "range1hi",
              "range2lo",
              "range2hi",
              "range3lo",
              "range3hi",
              "range4lo",
              "range4hi",
              "range5lo",
              "range5hi",
              "type")
colnames(parameter.space) <- varnames

# Export data to csv
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}

write.csv(parameter.space, "Parameter space from 2006HSE_25FEB2020.csv")

# Generate 'best estimates' csv input
best <- cbind(lnbmi.men$estimate[1],
          lnbmi.men$estimate[2],
          lnbmi.women$estimate[1],
          lnbmi.women$estimate[2],
          lnsbp.men$estimate[1],
          lnsbp.men$estimate[2],
          lnsbp.women$estimate[1],
          lnsbp.women$estimate[2],
          gammatcho.men$estimate[1],
          gammatcho.men$estimate[2],
          gammatcho.women$estimate[1],
          gammatcho.women$estimate[2],
          sum(hse3yr$x[hse3yr$sex==1&hse3yr$smok==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE),
          sum(hse3yr$x[hse3yr$sex==2&hse3yr$smok==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE),
          sum(hse3yr$x[hse3yr$sex==1&hse3yr$diab==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE),
          sum(hse3yr$x[hse3yr$sex==2&hse3yr$diab==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE),
          sum(hse3yr$x[hse3yr$sex==1&hse3yr$mi==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==1], na.rm = TRUE),
          sum(hse3yr$x[hse3yr$sex==2&hse3yr$mi==1], na.rm = TRUE)/sum(hse3yr$x[hse3yr$sex==2], na.rm = TRUE),
          age.bmi.mr.m.estimate,
          age.sbp.mr.m.estimate,
          age.tcho.mr.m.estimate,
          age.smok.mr.m.estimate,
          age.diab.mr.m.estimate,
          age.mi.mr.m.estimate,
          bmi.sbp.mr.m.estimate,
          bmi.tcho.mr.m.estimate,
          bmi.smok.mr.m.estimate,
          bmi.diab.mr.m.estimate,
          bmi.mi.mr.m.estimate,
          sbp.tcho.mr.m.estimate,
          sbp.smok.mr.m.estimate,
          sbp.diab.mr.m.estimate,
          sbp.mi.mr.m.estimate,
          tcho.smok.mr.m.estimate,
          tcho.diab.mr.m.estimate,
          tcho.mi.mr.m.estimate,
          smok.diab.mr.m.estimate,
          smok.mi.mr.m.estimate,
          diab.mi.mr.m.estimate,
          age.bmi.mr.f.estimate,
          age.sbp.mr.f.estimate,
          age.tcho.mr.f.estimate,
          age.smok.mr.f.estimate,
          age.diab.mr.f.estimate,
          age.mi.mr.f.estimate,
          bmi.sbp.mr.f.estimate,
          bmi.tcho.mr.f.estimate,
          bmi.smok.mr.f.estimate,
          bmi.diab.mr.f.estimate,
          bmi.mi.mr.f.estimate,
          sbp.tcho.mr.f.estimate,
          sbp.smok.mr.f.estimate,
          sbp.diab.mr.f.estimate,
          sbp.mi.mr.f.estimate,
          tcho.smok.mr.f.estimate,
          tcho.diab.mr.f.estimate,
          tcho.mi.mr.f.estimate,
          smok.diab.mr.f.estimate,
          smok.mi.mr.f.estimate,
          diab.mi.mr.f.estimate)

parameter.names <- c("bmi.pos.m",
                     "bmi.shp.m",
                     "bmi.pos.f",
                     "bmi.shp.f",
                     "sbp.pos.m",
                     "sbp.shp.m",
                     "sbp.pos.f",
                     "sbp.shp.f",
                     "tcho.pos.m",
                     "tcho.shp.m",
                     "tcho.pos.f",
                     "tcho.shp.f",
                     "smok.prev.m",
                     "smok.prev.f",
                     "diab.prev.m",
                     "diab.prev.f",
                     "mi.prev.m",
                     "mi.prev.f",
                     "age.bmi.mr.m",
                     "age.sbp.mr.m",
                     "age.tcho.mr.m",
                     "age.smok.mr.m",
                     "age.diab.mr.m",
                     "age.mi.mr.m",
                     "bmi.sbp.mr.m",
                     "bmi.tcho.mr.m",
                     "bmi.smok.mr.m",
                     "bmi.diab.mr.m",
                     "bmi.mi.mr.m",
                     "sbp.tcho.mr.m",
                     "sbp.smok.mr.m",
                     "sbp.diab.mr.m",
                     "sbp.mi.mr.m",
                     "tcho.smok.mr.m",
                     "tcho.diab.mr.m",
                     "tcho.mi.mr.m",
                     "smok.diab.mr.m",
                     "smok.mi.mr.m",
                     "diab.mi.mr.m",
                     "age.bmi.mr.f",
                     "age.sbp.mr.f",
                     "age.tcho.mr.f",
                     "age.smok.mr.f",
                     "age.diab.mr.f",
                     "age.mi.mr.f",
                     "bmi.sbp.mr.f",
                     "bmi.tcho.mr.f",
                     "bmi.smok.mr.f",
                     "bmi.diab.mr.f",
                     "bmi.mi.mr.f",
                     "sbp.tcho.mr.f",
                     "sbp.smok.mr.f",
                     "sbp.diab.mr.f",
                     "sbp.mi.mr.f",
                     "tcho.smok.mr.f",
                     "tcho.diab.mr.f",
                     "tcho.mi.mr.f",
                     "smok.diab.mr.f",
                     "smok.mi.mr.f",
                     "diab.mi.mr.f")
colnames(best) <- parameter.names

# Export data to csv
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}

write.csv(best, "Parameter space from 2006HSE_best estimates_25FEB2020.csv")

##################################################################################################
#                                                                                                #
# 4. Descriptions of trends in risk factors by age and year, stratified by sex                   #
#                                                                                                #
##################################################################################################

# Introduce variables for regression
hse.master$age2 <- hse.master$age^2
hse.master$year3 <- hse.master$year2 - 1992
hse.master$n <- 1

## BMI
# Eyeball the data, stratified by sex, changing by year and by age
# First, generate dataset of age-standardised mean of risk factors by year and sex
# NB: standardised to population in total hse.master dataset
awf <- function (sex, a1, a2) { # Function to produce age weights
  sum(na.omit(hse.master$n[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2]))/sum(na.omit(hse.master$n[hse.master$sex==sex]))
}
weight.m <- c(awf(1, 18, 25),
              awf(1, 25, 30),
              awf(1, 30, 35),
              awf(1, 35, 40),
              awf(1, 40, 45),
              awf(1, 45, 50),
              awf(1, 50, 55),
              awf(1, 55, 60),
              awf(1, 60, 65),
              awf(1, 65, 70),
              awf(1, 70, 75),
              awf(1, 75, 80),
              awf(1, 80, 85),
              awf(1, 85, 150))
weight.f <- c(awf(2, 18, 25),
              awf(2, 25, 30),
              awf(2, 30, 35),
              awf(2, 35, 40),
              awf(2, 40, 45),
              awf(2, 45, 50),
              awf(2, 50, 55),
              awf(2, 55, 60),
              awf(2, 60, 65),
              awf(2, 65, 70),
              awf(2, 70, 75),
              awf(2, 75, 80),
              awf(2, 80, 85),
              awf(2, 85, 150))

mwf.bmi <- function (sex, a1, a2, year) { # Function to produce mean of bmi
  mean(na.omit(hse.master$bmi[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2&hse.master$year2==year]))
}
mwf.sbp <- function (sex, a1, a2, year) { # Function to produce mean of sbp
  mean(na.omit(hse.master$sbp[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2&hse.master$year2==year]))
}
mwf.chol <- function (sex, a1, a2, year) { # Function to produce mean of chol
  mean(na.omit(hse.master$chol[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2&hse.master$year2==year]))
}
mwf.diab <- function (sex, a1, a2, year) { # Function to produce mean of diabetes
  mean(na.omit(hse.master$diabetes[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2&hse.master$year2==year]))
}
mwf.smok <- function (sex, a1, a2, year) { # Function to produce mean of smok
  mean(na.omit(hse.master$smok[hse.master$sex==sex&hse.master$age>=a1&hse.master$age<a2&hse.master$year2==year]))
}

for (i in 1991:2014){
  assign(paste0("bmi.m.",i),sum(c(mwf.bmi(1, 18, 25, i),
                                  mwf.bmi(1, 25, 30, i),
                                  mwf.bmi(1, 30, 35, i),
                                  mwf.bmi(1, 35, 40, i),
                                  mwf.bmi(1, 40, 45, i),
                                  mwf.bmi(1, 45, 50, i),
                                  mwf.bmi(1, 50, 55, i),
                                  mwf.bmi(1, 55, 60, i),
                                  mwf.bmi(1, 60, 65, i),
                                  mwf.bmi(1, 65, 70, i),
                                  mwf.bmi(1, 70, 75, i),
                                  mwf.bmi(1, 75, 80, i),
                                  mwf.bmi(1, 80, 85, i),
                                  mwf.bmi(1, 85, 150, i))*weight.m))
  assign(paste0("bmi.f.",i),sum(c(mwf.bmi(2, 18, 25, i),
                                  mwf.bmi(2, 25, 30, i),
                                  mwf.bmi(2, 30, 35, i),
                                  mwf.bmi(2, 35, 40, i),
                                  mwf.bmi(2, 40, 45, i),
                                  mwf.bmi(2, 45, 50, i),
                                  mwf.bmi(2, 50, 55, i),
                                  mwf.bmi(2, 55, 60, i),
                                  mwf.bmi(2, 60, 65, i),
                                  mwf.bmi(2, 65, 70, i),
                                  mwf.bmi(2, 70, 75, i),
                                  mwf.bmi(2, 75, 80, i),
                                  mwf.bmi(2, 80, 85, i),
                                  mwf.bmi(2, 85, 150, i))*weight.f))
  assign(paste0("sbp.m.",i),sum(c(mwf.sbp(1, 18, 25, i),
                                  mwf.sbp(1, 25, 30, i),
                                  mwf.sbp(1, 30, 35, i),
                                  mwf.sbp(1, 35, 40, i),
                                  mwf.sbp(1, 40, 45, i),
                                  mwf.sbp(1, 45, 50, i),
                                  mwf.sbp(1, 50, 55, i),
                                  mwf.sbp(1, 55, 60, i),
                                  mwf.sbp(1, 60, 65, i),
                                  mwf.sbp(1, 65, 70, i),
                                  mwf.sbp(1, 70, 75, i),
                                  mwf.sbp(1, 75, 80, i),
                                  mwf.sbp(1, 80, 85, i),
                                  mwf.sbp(1, 85, 150, i))*weight.m))
  assign(paste0("sbp.f.",i),sum(c(mwf.sbp(2, 18, 25, i),
                                  mwf.sbp(2, 25, 30, i),
                                  mwf.sbp(2, 30, 35, i),
                                  mwf.sbp(2, 35, 40, i),
                                  mwf.sbp(2, 40, 45, i),
                                  mwf.sbp(2, 45, 50, i),
                                  mwf.sbp(2, 50, 55, i),
                                  mwf.sbp(2, 55, 60, i),
                                  mwf.sbp(2, 60, 65, i),
                                  mwf.sbp(2, 65, 70, i),
                                  mwf.sbp(2, 70, 75, i),
                                  mwf.sbp(2, 75, 80, i),
                                  mwf.sbp(2, 80, 85, i),
                                  mwf.sbp(2, 85, 150, i))*weight.f))
  assign(paste0("chol.m.",i),sum(c(mwf.chol(1, 18, 25, i),
                                  mwf.chol(1, 25, 30, i),
                                  mwf.chol(1, 30, 35, i),
                                  mwf.chol(1, 35, 40, i),
                                  mwf.chol(1, 40, 45, i),
                                  mwf.chol(1, 45, 50, i),
                                  mwf.chol(1, 50, 55, i),
                                  mwf.chol(1, 55, 60, i),
                                  mwf.chol(1, 60, 65, i),
                                  mwf.chol(1, 65, 70, i),
                                  mwf.chol(1, 70, 75, i),
                                  mwf.chol(1, 75, 80, i),
                                  mwf.chol(1, 80, 85, i),
                                  mwf.chol(1, 85, 150, i))*weight.m))
  assign(paste0("chol.f.",i),sum(c(mwf.chol(2, 18, 25, i),
                                  mwf.chol(2, 25, 30, i),
                                  mwf.chol(2, 30, 35, i),
                                  mwf.chol(2, 35, 40, i),
                                  mwf.chol(2, 40, 45, i),
                                  mwf.chol(2, 45, 50, i),
                                  mwf.chol(2, 50, 55, i),
                                  mwf.chol(2, 55, 60, i),
                                  mwf.chol(2, 60, 65, i),
                                  mwf.chol(2, 65, 70, i),
                                  mwf.chol(2, 70, 75, i),
                                  mwf.chol(2, 75, 80, i),
                                  mwf.chol(2, 80, 85, i),
                                  mwf.chol(2, 85, 150, i))*weight.f))
  assign(paste0("diab.m.",i),sum(c(mwf.diab(1, 18, 25, i),
                                  mwf.diab(1, 25, 30, i),
                                  mwf.diab(1, 30, 35, i),
                                  mwf.diab(1, 35, 40, i),
                                  mwf.diab(1, 40, 45, i),
                                  mwf.diab(1, 45, 50, i),
                                  mwf.diab(1, 50, 55, i),
                                  mwf.diab(1, 55, 60, i),
                                  mwf.diab(1, 60, 65, i),
                                  mwf.diab(1, 65, 70, i),
                                  mwf.diab(1, 70, 75, i),
                                  mwf.diab(1, 75, 80, i),
                                  mwf.diab(1, 80, 85, i),
                                  mwf.diab(1, 85, 150, i))*weight.m))
  assign(paste0("diab.f.",i),sum(c(mwf.diab(2, 18, 25, i),
                                  mwf.diab(2, 25, 30, i),
                                  mwf.diab(2, 30, 35, i),
                                  mwf.diab(2, 35, 40, i),
                                  mwf.diab(2, 40, 45, i),
                                  mwf.diab(2, 45, 50, i),
                                  mwf.diab(2, 50, 55, i),
                                  mwf.diab(2, 55, 60, i),
                                  mwf.diab(2, 60, 65, i),
                                  mwf.diab(2, 65, 70, i),
                                  mwf.diab(2, 70, 75, i),
                                  mwf.diab(2, 75, 80, i),
                                  mwf.diab(2, 80, 85, i),
                                  mwf.diab(2, 85, 150, i))*weight.f))
  assign(paste0("smok.m.",i),sum(c(mwf.smok(1, 18, 25, i),
                                  mwf.smok(1, 25, 30, i),
                                  mwf.smok(1, 30, 35, i),
                                  mwf.smok(1, 35, 40, i),
                                  mwf.smok(1, 40, 45, i),
                                  mwf.smok(1, 45, 50, i),
                                  mwf.smok(1, 50, 55, i),
                                  mwf.smok(1, 55, 60, i),
                                  mwf.smok(1, 60, 65, i),
                                  mwf.smok(1, 65, 70, i),
                                  mwf.smok(1, 70, 75, i),
                                  mwf.smok(1, 75, 80, i),
                                  mwf.smok(1, 80, 85, i),
                                  mwf.smok(1, 85, 150, i))*weight.m))
  assign(paste0("smok.f.",i),sum(c(mwf.smok(2, 18, 25, i),
                                  mwf.smok(2, 25, 30, i),
                                  mwf.smok(2, 30, 35, i),
                                  mwf.smok(2, 35, 40, i),
                                  mwf.smok(2, 40, 45, i),
                                  mwf.smok(2, 45, 50, i),
                                  mwf.smok(2, 50, 55, i),
                                  mwf.smok(2, 55, 60, i),
                                  mwf.smok(2, 60, 65, i),
                                  mwf.smok(2, 65, 70, i),
                                  mwf.smok(2, 70, 75, i),
                                  mwf.smok(2, 75, 80, i),
                                  mwf.smok(2, 80, 85, i),
                                  mwf.smok(2, 85, 150, i))*weight.f))
}

# Create year collapsed dataset
hse.year.collapse <- NULL
for (i in 1991:2014) {
  Z <- c(i, 
         get(paste0("bmi.m.",i)),
         get(paste0("bmi.f.",i)),
         get(paste0("sbp.m.",i)),
         get(paste0("sbp.f.",i)),
         get(paste0("chol.m.",i)),
         get(paste0("chol.f.",i)),
         get(paste0("diab.m.",i)),
         get(paste0("diab.f.",i)),
         get(paste0("smok.m.",i)),
         get(paste0("smok.f.",i)))
  hse.year.collapse <- rbind(hse.year.collapse,Z)
}

names <- c("year",
           "bmi.m",
           "bmi.f",
           "sbp.m",
           "sbp.f",
           "chol.m",
           "chol.f",
           "diab.m",
           "diab.f",
           "smok.m",
           "smok.f")

colnames(hse.year.collapse) <- names

hse.year.collapse <- as.data.frame(hse.year.collapse)

# Now produce dataset collapsed on five-year age group
# Produce a five year age variable in hse.master
agefunction <- function (age){
  if (is.na(age)==TRUE) { NA }
  else if (age<25)  { 1 }
  else if (age<85)  { floor((age-15)/5) }
  else              { 14 }
}

hse.master$agegroup <- as.numeric(lapply(hse.master$age, agefunction))

# Collapse on agegroup

hse.age.collapse.m <- ddply(hse.master[hse.master$sex==1,], 
                            .(agegroup), summarize,
                            bmi=mean(bmi, na.rm = TRUE),
                            sbp=mean(sbp, na.rm = TRUE),
                            chol=mean(chol, na.rm=TRUE),
                            diab=mean(diabetes, na.rm=TRUE),
                            smok=mean(smok, na.rm=TRUE))

hse.age.collapse.f <- ddply(hse.master[hse.master$sex==2,], 
                            .(agegroup), summarize,
                            bmi=mean(bmi, na.rm = TRUE),
                            sbp=mean(sbp, na.rm = TRUE),
                            chol=mean(chol, na.rm=TRUE),
                            diab=mean(diabetes, na.rm=TRUE),
                            smok=mean(smok, na.rm=TRUE))

# Plot relationship between variables and age and year

### BMI, MEN

# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$bmi) # Definite quadratic fit
plot (hse.year.collapse$year,
      hse.year.collapse$bmi.m) # Maybe linear, but probably quadratic

# Generate year variable with zero at 1991
hse.master$year3 <- as.numeric(hse.master$year2) - 1991

# Generate quadratic variables
hse.master$agesq <- hse.master$age^2
hse.master$yearsq <- hse.master$year3^2

# Model 1: quadratic age, linear year
m1 <- lm(bmi ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==1,])
# Model 2: quadratic age, quadratic year
m2 <- lm(bmi ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==1,])
# LR test to see if better fit
lrtest(m2, m1) # p<0.0001

# Check whether older age groups have significantly different trends over time than younger
# First introduce binary 'old' variable, where old is defined as 60 and over.
oldfunction <- function(x){
  if (is.na(x)==TRUE) { NA }
  else if (x>=60)     { 1 }
  else                { 0 }
}
hse.master$old <- lapply(hse.master$age, oldfunction)
# Now generate interaction terms
hse.master$oldyearint   <- as.numeric(hse.master$old)*as.numeric(hse.master$year3)
hse.master$oldyearsqint <- as.numeric(hse.master$old)*as.numeric(hse.master$yearsq)
# New model with interaction terms
m3 <- lm(bmi ~ age + agesq + year3 + yearsq + oldyearint + oldyearsqint,
         data = hse.master[hse.master$sex==1,])
# LR test to see if better fit
lrtest(m3, m2) # p<0.0001

# Try a model with a continuously changing interaction term
# First generate interaction terms
hse.master$ageyearint   <- as.numeric(hse.master$age)*as.numeric(hse.master$year3)
hse.master$ageyearsqint <- as.numeric(hse.master$age)*as.numeric(hse.master$yearsq)
# Then model and LR test
m4 <- lm(bmi ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==1,])
lrtest(m4, m2) # p<0.0001
# OK - but the interaction terms are not significant.

# Try a linear model with interaction terms
m5 <- lm(bmi ~ age + agesq + year3 + ageyearint,
         data = hse.master[hse.master$sex==1,])
lrtest(m5, m1) # p<0.0001

### BMI, WOMEN

# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.f$agegroup, 
      hse.age.collapse.f$bmi) # Definite quadratic fit, with some funky data at old ages
plot (hse.year.collapse$year,
      hse.year.collapse$bmi.f) # Maybe linear, but probably quadratic

# Try different models and check for differences with LRtests
f1 <- lm(bmi ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==2,])
f2 <- lm(bmi ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==2,])
lrtest(f2, f1) # p<0.0001
f3 <- lm(bmi ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==2,])
lrtest(f3, f2) # p = 0.1849

### SBP, MEN

# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$sbp) # Definite quadratic fit, with near flat trend at old age
plot (hse.year.collapse$year,
      hse.year.collapse$sbp.m) # Maybe linear, but probably quadratic

# Try different models and check for differences with LRtests
m1 <- lm(sbp ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==1,])
m2 <- lm(sbp ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==1,]) # Note that the quadratic term is +ve, contrary to the plots
                                                # Model is fitting to the young data
lrtest(m2, m1) # p<0.0001
m3 <- lm(sbp ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==1,])
lrtest(m3, m2) # p<0.0001

# Need to rerun the quadratic model with the young ages removed. This will force the quadratic curve
# to go through the data at the old age points, and we can force the blood pressure to stay
# constant at early ages.
m1 <- lm(sbp ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==1&hse.master$age>=40,])
m2 <- lm(sbp ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==1&hse.master$age>=40,]) 
lrtest(m2, m1) # p<0.0001
# Test for interactions
m3 <- lm(sbp ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==1&hse.master$age>=40,])
lrtest(m3, m2) # p<0.0001

# SBP WOMEN
# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.f$agegroup, 
      hse.age.collapse.f$sbp) # Definite quadratic fit, with near flat trend at old age
plot (hse.year.collapse$year,
      hse.year.collapse$sbp.f) # Maybe linear, but probably quadratic. Very similar to men. Take the same analysis approach.

f1 <- lm(sbp ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==2&hse.master$age>=40,])
f2 <- lm(sbp ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==2&hse.master$age>=40,]) 
lrtest(f2, f1) # p<0.0001

# TCHO MEN
# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$chol) # Definite quadratic fit
plot (hse.year.collapse$year,
      hse.year.collapse$chol.m) # Patchy, but probably linear trend down

m1 <- lm(chol ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==1,])
m2 <- lm(chol ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==1,]) 
lrtest(m2, m1) # p=0.2466 indicating no evidence of quadratic trend in year

# Try interaction with age
m3 <- lm(chol ~ age + agesq + year3 + ageyearint,
         data = hse.master[hse.master$sex==1,])
lrtest(m3, m1) # p<0.001

# The interaction one looks ok, but some projections go below zero. Try an interaction model with a quadratic fit for time
m4 <- lm(chol ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==1,])
lrtest(m4, m3) # p=0.00357. But the three of the exposures are barely significant.
# Go for the version without interactions, as with interactions leads to implausible projections.

# TCHO WOMEN
# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.f$agegroup, 
      hse.age.collapse.f$chol) # Definite quadratic fit
plot (hse.year.collapse$year,
      hse.year.collapse$chol.f) # Patchy, but probably linear trend down

f1 <- lm(chol ~ age + agesq + year3, 
         data = hse.master[hse.master$sex==2,])
f2 <- lm(chol ~ age + agesq + year3 + yearsq, 
         data = hse.master[hse.master$sex==2,]) 
lrtest(f2, f1) # p<0.001

# Try interactions
f3 <- lm(chol ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
         data = hse.master[hse.master$sex==2,])
lrtest(f3, f2) # p<0.001

# SMOKING MEN
# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$smok) # Looks like a linear decline
plot (hse.year.collapse$year,
      hse.year.collapse$smok.m) # Looks like a linear decline

m1 <- glm(smok ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
m2 <- glm(smok ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m2, m1) # p<0.001

m3 <- glm(smok ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m3, m2) # p = 0.5558
m4 <- glm(smok ~ age + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m4, m1) # p = 0.7225 No evidence of quadratic trend for time

# Try an interacton term
m5 <- glm(smok ~ age + agesq + year3 + ageyearint,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m5, m2) # p<0.001, but year coefficient no longer significant

# SMOKING WOMEN
# Observe the data 
par(mfrow=c(2,1))
plot (hse.age.collapse.f$agegroup, 
      hse.age.collapse.f$smok) # Looks like a linear decline
plot (hse.year.collapse$year,
      hse.year.collapse$smok.f) # Looks like a linear decline

f1 <- glm(smok ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==2,])
f2 <- glm(smok ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==2,])
lrtest(f2, f1) # p<0.001

f3 <- glm(smok ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==2,])
lrtest(f3, f2) # p = 0.00222
f4 <- glm(smok ~ age + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==2,])
lrtest(f4, f1) # p = 0.00579. Some evidence of quadratic trend for time. Model f3 is best.

# Try an interacton term
f5 <- glm(smok ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
          family = "binomial",
          data = hse.master[hse.master$sex==2,])
lrtest(f5, f3) # p<0.001, but year coefficient no longer significant


# DIABETES MEN
# Observe the data
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$diab) # Looks like a quadatric increase. May need to cut off young data points
plot (hse.year.collapse$year,
      hse.year.collapse$diab.m) # Looks like a linear decline

m1 <- glm(diabetes ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
m2 <- glm(diabetes ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m2, m1) # p<0.001
m3 <- glm(diabetes ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==1,])
lrtest(m3, m2) # p<0.001
# Try models from 1994 onwards, since weird spike at 1993
par(mfrow=c(2,1))
plot (hse.age.collapse.m$agegroup, 
      hse.age.collapse.m$diab) # Looks like a quadatric increase. May need to cut off young data points
plot (hse.year.collapse$year[hse.year.collapse$year>=1994],
      hse.year.collapse$diab.m[hse.year.collapse$year>=1994]) # Looks like a linear decline
m4 <- glm(diabetes ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$year3>=3,])
m5 <- glm(diabetes ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$year3>=3,])
lrtest(m5, m4) # p<0.001
m6 <- glm(diabetes ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$year3>=3,])
lrtest(m6, m5) # p=0.001606. This is better
# The quadratic for year is positive, suggesting an accelerating trend, which is not observed in the data,
# unless the outliers have a massive impact. So drop the outliers.
diaboutlier <- function (x) {
  if (x==1993|x==1999|x==2004|x==2005) { 1 }
  else { 0 }
}
hse.master$diabout <- as.numeric(lapply(hse.master$year2, diaboutlier))
m7 <- glm(diabetes ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$diabout!=1,])
m8 <- glm(diabetes ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$diabout!=1,])
lrtest(m8, m7) # p<0.001
m9 <- glm(diabetes ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==1&hse.master$diabout!=1,])
lrtest(m9, m8) # p = 0.02282

# Try an interaction term
m10 <- glm(diabetes ~ age + agesq + year3 + yearsq + ageyearint + ageyearsqint,
           family = "binomial",
           data = hse.master[hse.master$sex==1&hse.master$diabout!=1,])
lrtest(m10, m9) # p = 0.1726. No evidence of improved fit.

# DIABETES WOMEN
# Observe the data
par(mfrow=c(2,1))
plot (hse.age.collapse.f$agegroup, 
      hse.age.collapse.f$diab) # Looks like a quadatric increase. May need to cut off young data points
plot (hse.year.collapse$year[hse.year.collapse$year>=1994],
      hse.year.collapse$diab.f[hse.year.collapse$year>=1994]) # Need to remove outliers

f1 <- glm(diabetes ~ age + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==2&hse.master$diabout!=1,])
f2 <- glm(diabetes ~ age + agesq + year3,
          family = "binomial",
          data = hse.master[hse.master$sex==2&hse.master$diabout!=1,])
lrtest(f2, f1) ~ p<0.001

f3 <- glm(diabetes ~ age + agesq + year3 + yearsq,
          family = "binomial",
          data = hse.master[hse.master$sex==2&hse.master$diabout!=1,])
lrtest(f3, f2) # p = 0.688 - No evidence of quadratic trend in year

# Try for an interaction term
f4 <- glm(diabetes ~ age + agesq + year3 + ageyearint,
          family = "binomial",
          data = hse.master[hse.master$sex==2&hse.master$diabout!=1,])
lrtest(f4, f2) # p = 0.3956. No evidence of interaction term.

# MI - men
f5 <- glm(mi ~ age + bmi + sbp + chol + smok + diabetes,
          family = "binomial",
          data = hse3yr[hse3yr$sex==1,])

hse3yr$agesq <- hse3yr$age*hse3yr$age

hse3yr$agecub <- hse3yr$agesq*hse3yr$age

f6 <- glm(mi ~ age + agesq,# + bmi + sbp + chol + smok + diabetes,
          family = "binomial",
          data = hse3yr[hse3yr$sex==1,])
lrtest(f6,f5) # Significantly better fit.

f7 <- glm(mi ~ age + agesq + agecub + bmi + sbp + chol + smok + diabetes,
          family = "binomial",
          data = hse3yr[hse3yr$sex==1,])
lrtest(f7,f6) # Not a better fit.

f8 <- glm(mi ~ age + agesq,# + bmi + sbp + chol + smok + diabetes,
          family = "binomial",
          data = hse3yr[hse3yr$sex==2,])

# Save f6 and f8 for applying in agent fill
coefficients <- cbind(f6$coefficients,f8$coefficients)

# Set directory for export
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}

write.csv(coefficients, "miprev_coefficients_28NOV2019.csv") # NB: Needs updating each time
