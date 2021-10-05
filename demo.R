### demo ###

# Last updated: 12th October 2018

# There are three aims of this R script:
# 1. To produce estimates of 'other mortality' rates, stratified by sex and smoothly varying by age and calendar year
# 2. To produce estimates of migration rates, stratified by sex and varying by age and calendar year
# 3. To produce estimates of fertility rates by calendar year

# Note that this R script relies on running the 'demo_preparing dataset.R' script to produce the initial dataset. See note
# below when the fert.mort dataset is read.

library("plyr")

work <- "desktop" # Choice between 'desktop' and 'laptop'. Ensures script uses correct directories.

# Load dataset
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Demographic trends")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Demographic trends")
  }  else { "Error" }
}

fert.mort <- read.csv("pop.mort04AUG2017.csv") # To update the input file, need to run 'demo_preparing dataset.R'


##################################################################################################
#                                                                                                #
# 1. Estimating and projecting 'other mortality' rates                                           #
#                                                                                                #
##################################################################################################

# New variable that is simply total mortality rate minus MI mortality rate
fert.mort$oth.mort.p <- fert.mort$mort.p - fert.mort$mimort.p

# Check plots of age against other mortality rate (as percentage of population) in 1992 and 2015
age.r <- c(2.5,        # This is the age variable for regression analyses and plots. Under 1s ignored, since we already have single
           (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
           85)         # For the plots below, this is arbitrarily set as 85.

par(mfrow=c(2,2))
plot(age.r, 
     fert.mort$oth.mort.p[(fert.mort$X%%5==3&fert.mort$year==1992&fert.mort$sex==1)|  # %% is modulo, takes every third reading
                            (fert.mort$age==85&fert.mort$year==1992&fert.mort$sex==1)], # This gets the final estimate for over 85s
     xlab = "Age", ylab = "Other mort rate, 1992, men")
plot(age.r, 
     fert.mort$oth.mort.p[(fert.mort$X%%5==3&fert.mort$year==1992&fert.mort$sex==2)|  
                            (fert.mort$age==85&fert.mort$year==1992&fert.mort$sex==2)],
     xlab = "Age", ylab = "Other mort rate, 1992, women") 
plot(age.r, 
     fert.mort$oth.mort.p[(fert.mort$X%%5==3&fert.mort$year==2015&fert.mort$sex==1)|  
                            (fert.mort$age==85&fert.mort$year==2015&fert.mort$sex==1)],
     xlab = "Age", ylab = "Other mort rate, 2015, men") 
plot(age.r, 
     fert.mort$oth.mort.p[(fert.mort$X%%5==3&fert.mort$year==2015&fert.mort$sex==2)|  
                            (fert.mort$age==85&fert.mort$year==2015&fert.mort$sex==2)], 
     xlab = "Age", ylab = "Other mort rate, 2015, women") 

# Looks exponential in shape, so can run exponential fit. Not clear what average age should be taken for the oldest group. Code
# below estimates the best choice to fit the data

# To find the parameter (y) for each year which minimises the square difference between predicted and measured rates:
sum.diff.men.by.j   <- c() # To store the results of the big loop initiated below
sum.diff.women.by.j <- c()

for (j in 65:95) {
  
  y <- j                 # This sets the age for oldest age group
  
  age.r <- c(2.5,        # This will be the age variable for regression analyses. Under 1s ignored, since we already have single
             (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
             y)          # See above.
  
  
  
  # Run exponential fitting models
  fert.mort$lnoth.mort.p                    <- log(fert.mort$oth.mort.p) # Produce log-transformed variable
  
  fert.mort$age.for.pred                    <- fert.mort$age  # This introduces a new variable with age 85 replaced with y as 
  fert.mort$age.for.pred[fert.mort$age==85] <- y              # average age for the over 85 group for the predicted rates.
  
  fert.mort$pred.oth.mort                   <- c(rep(0,4128)) # Dummy variable to be replaced with predicted rates in the loop below.
  fert.mort$cons                            <- c(rep(1,4128)) # Column of 1s for the prediction process.
  
  
  # Create reduced dataset for the regressions
  regress.db  <- fert.mort[fert.mort$age%%5==3|fert.mort$age.for.pred==y,] # NB: last comma is to show I want all variables
  regressions <- dlply(regress.db,                        # This runs a series of functions across diffent elements of regress.db
                       .(year, sex),                      # Split by year and sex
                       lm,                                # Linear models (i.e. regressions) are run
                       formula = lnoth.mort.p ~ age.for.pred) # To this formula
  coefs       <- ldply(regressions, coef)                 # This captures the results of the 48 regressions in a dataset
  
  for (i in 1992:2015) {
    fert.mort$pred.oth.mort[fert.mort$year==i&fert.mort$sex==1] <-              # Replacing pred.mort by year and sex group
      exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==1],          # Make a 86x2 matrix of 1s and ages
                   fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==1]), # 'exp' is to unlog model predictions
                 ncol = 2, nrow = 86)
          %*%matrix(c(coefs[(2*(i-1992)+1),3],coefs[(2*(i-1992)+1),4]),         # Multiply by 2x1 matrix of year and sex
                    ncol = 1, nrow = 2))                                        # specific regression coefficients.
    
    fert.mort$pred.oth.mort[fert.mort$year==i&fert.mort$sex==2] <-
      exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==2],
                   fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==2]),
                 ncol = 2, nrow = 86)
          %*%matrix(c(coefs[(2*(i-1992)+2),3],coefs[(2*(i-1992)+2),4]),
                    ncol = 1, nrow = 2))
  }
  
  # Aim is to select y to minimise the difference between predicted and measured mortality rates up to 83, and then use that.
  sum.diff.men        <- c() # To store the results of the small loop below
  sum.diff.women      <- c()
  for (i in 1992:2015){ # This is to calculate the sum of the square of differences between prediction and measure
    z1 <- sum(fert.mort$pred.oth.mort[fert.mort$age%%5==3&fert.mort$sex==1&fert.mort$year==i]-
                fert.mort$oth.mort.p[fert.mort$age%%5==3&fert.mort$sex==1&fert.mort$year==i])^2
    sum.diff.men <- append(sum.diff.men, z1)
    z2 <- sum(fert.mort$pred.oth.mort[fert.mort$age%%5==3&fert.mort$sex==2&fert.mort$year==i]-
                fert.mort$oth.mort.p[fert.mort$age%%5==3&fert.mort$sex==2&fert.mort$year==i])^2
    sum.diff.women <- append(sum.diff.women, z2)
  }
  
  sum.diff.men.by.j   <- rbind(sum.diff.men.by.j,   sum.diff.men)
  sum.diff.women.by.j <- rbind(sum.diff.women.by.j, sum.diff.women)
  
}

sum.diff.men.by.j   <- as.data.frame(sum.diff.men.by.j)
sum.diff.women.by.j <- as.data.frame(sum.diff.women.by.j)

# Find minimum sum of squared difference for each year
min.men   <- c()
min.women <- c()
for (i in 1:24){
  z1 <- which.min(sum.diff.men.by.j[,i])+64
  min.men <- append(min.men, z1)
  z2 <- which.min(sum.diff.women.by.j[,i])+64
  min.women <- append(min.women, z2)
}

# Plot the results that show where the best exponential fit is
year <- 1992:2015
par(mfrow=c(1,2))
plot(year, min.men,   ylab = "Minimum squared diff by j, men")
plot(year, min.women, ylab = "Minimum squared diff by j, women")
# This demonstrates that for both men and women, the old age group should not be more than 85, but for many years you get a 
# better fit of the model if the oldest age group is less than 85. But that doesn't make sense, so from herein the oldest age is
# set as 85.

## WARNING: THE PLOTS SUGGEST THAT THE VALUE FOR Y INCREASES AS TIME INCREASES (WHICH MAKES SENSE), WHICH MAY MEAN THAT USING
## 85 AS A VALUE FOR Y WILL BE WRONG FOR PROJECTIONS.

# Set the predicted other mortality rates with oldest age group set as 85
y <- 85
age.r <- c(2.5,        # This will be the age variable for regression analyses. Under 1s ignored, since we already have single
           (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
           y)          # See above.
fert.mort$age.for.pred[fert.mort$age==85] <- y              # average age for the over 85 group for the predicted rates.
fert.mort$pred.oth.mort                   <- c(rep(0,4128)) # Dummy variable to be replaced with predicted rates in the loop below.
fert.mort$cons                            <- c(rep(1,4128)) # Column of 1s for the prediction process.
# Create reduced dataset for the regressions
regress.db <- fert.mort[fert.mort$age%%5==3|fert.mort$age.for.pred==y,] # NB: last comma is to show I want all variables
regressions <- dlply(regress.db,                        # This runs a series of functions across diffent elements of regress.db
                     .(year, sex),                      # Split by year and sex
                     lm,                                # Linear models (i.e. regressions) are run
                     formula = lnoth.mort.p ~ age.for.pred) # To this formula
coefs       <- ldply(regressions, coef)                 # This captures the results of the 48 regressions in a dataset

for (i in 1992:2015) {
  fert.mort$pred.oth.mort[fert.mort$year==i&fert.mort$sex==1] <-              # Replacing pred.mort by year and sex group
    exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==1],          # Make a 86x2 matrix of 1s and ages
                 fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==1]), # 'exp' is to unlog model predictions
               ncol = 2, nrow = 86)
        %*%matrix(c(coefs[(2*(i-1992)+1),3],coefs[(2*(i-1992)+1),4]),         # Multiply by 2x1 matrix of year and sex
                  ncol = 1, nrow = 2))                                        # specific regression coefficients.
  
  fert.mort$pred.oth.mort[fert.mort$year==i&fert.mort$sex==2] <-
    exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==2],
                 fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==2]),
               ncol = 2, nrow = 86)
        %*%matrix(c(coefs[(2*(i-1992)+2),3],coefs[(2*(i-1992)+2),4]),
                  ncol = 1, nrow = 2))
}

# Check predicted results against measured results for men and women in 1992 and 2015
# Red = measured results; Green = predicted results
par(mfrow=c(2,2))
plot(regress.db$age.for.pred[regress.db$year==1992&regress.db$sex==1], # Men, 1992
     regress.db$oth.mort.p[regress.db$year==1992&regress.db$sex==1], 
     col = "red", xlab = "Age", ylab = "Other mort rate, 1992, men")
points(fert.mort$age.for.pred[fert.mort$year==1992&fert.mort$sex==1],
       fert.mort$pred.oth.mort[fert.mort$year==1992&fert.mort$sex==1],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==1992&regress.db$sex==2], # Women, 1992
     regress.db$oth.mort.p[regress.db$year==1992&regress.db$sex==2], 
     col = "red", xlab = "Age", ylab = "Other mort rate, 1992, women")
points(fert.mort$age.for.pred[fert.mort$year==1992&fert.mort$sex==2],
       fert.mort$pred.oth.mort[fert.mort$year==1992&fert.mort$sex==2],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==2015&regress.db$sex==1], # Men, 2015
     regress.db$oth.mort.p[regress.db$year==2015&regress.db$sex==1], 
     col = "red", xlab = "Age", ylab = "Other mort rate, 2015, men")
points(fert.mort$age.for.pred[fert.mort$year==2015&fert.mort$sex==1],
       fert.mort$pred.oth.mort[fert.mort$year==2015&fert.mort$sex==1],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==2015&regress.db$sex==2], # Women, 2015
     regress.db$oth.mort.p[regress.db$year==2015&regress.db$sex==2], 
     col = "red", xlab = "Age", ylab = "Other mort rate, 2015, women")
points(fert.mort$age.for.pred[fert.mort$year==2015&fert.mort$sex==2],
       fert.mort$pred.oth.mort[fert.mort$year==2015&fert.mort$sex==2],
       col = "green")
# All close up to 83 but suggest that the oldest age mortality rate is not reached until much after 85 (especially in later years).
# For 2015 women, the rate reaches 0.15 at age 99.

# Plot of predicted other mortality rate by age, year and sex, extended to age 100
X                      <- 1:4848
modelled.mort          <- data.frame(X)                    # Place holder to set up new dataframe
modelled.mort$age      <- rep(0:100,48)                    # Introducing the age,
modelled.mort$sex      <- rep(c(rep(1,101),rep(2,101)),24) # sex,
modelled.mort$year     <- floor((X-1)/202)+1992            # and year variables.
modelled.mort$cons     <- 1                                # And a row of ones for applying the model (see below)
modelled.mort$oth.mort <- 1                                # Place holder replaced below

for (i in 1992:2015) {
  modelled.mort$oth.mort[modelled.mort$year==i&modelled.mort$sex==1] <-          # Generating other mortality by year and sex group
    exp(matrix(c(modelled.mort$cons[modelled.mort$year==i&modelled.mort$sex==1], # Make a 101x2 matrix of 1s and ages
                 modelled.mort$age[modelled.mort$year==i&modelled.mort$sex==1]), # 'exp' is to unlog model predictions
               ncol = 2, nrow = 101)
        %*%matrix(c(coefs[(2*(i-1992)+1),3],coefs[(2*(i-1992)+1),4]),         # Multiply by 2x1 matrix of year and sex
                  ncol = 1, nrow = 2))                                        # specific regression coefficients.
  
  modelled.mort$oth.mort[modelled.mort$year==i&modelled.mort$sex==2] <-
    exp(matrix(c(modelled.mort$cons[modelled.mort$year==i&modelled.mort$sex==2],
                 modelled.mort$age[modelled.mort$year==i&modelled.mort$sex==2]),
               ncol = 2, nrow = 101)
        %*%matrix(c(coefs[(2*(i-1992)+2),3],coefs[(2*(i-1992)+2),4]),
                  ncol = 1, nrow = 2))
}

# Plots
par(mfrow=c(1,2))
plot(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1992],
     modelled.mort$oth.mort[modelled.mort$sex==1&modelled.mort$year==1992],
     xlab = "Age", ylab = "Other mortality rate, men", col = 92)
for (i in 1992:2015){
  points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==i],
         modelled.mort$oth.mort[modelled.mort$sex==1&modelled.mort$year==i],
         col = i-1900)
  }
plot(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1992],
     modelled.mort$oth.mort[modelled.mort$sex==2&modelled.mort$year==1992],
     xlab = "Age", ylab = "Other mortality rate, women", col = 92)
for (i in 1992:2015){
  points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==i],
         modelled.mort$oth.mort[modelled.mort$sex==2&modelled.mort$year==i],
         col = i-1900)
}

# The plots suggest that the 'other mortality rates' are reducing linearly by age (within the exponential increase by age)

# Try the following model: Other mortality rate = exp(a + b1*age + b2*year)
# Fit using the fert.mort dataset and compare using the modelled.mort dataset
# Fitting the relationship
regress.db   <- fert.mort[fert.mort$age%%5==3|fert.mort$age.for.pred==y,] # NB: last comma is to show I want all variables
regressions2 <- dlply(regress.db,                        # This runs a series of functions across diffent elements of regress.db
                     .(sex),                            # Split by sex
                     lm,                                # Linear models (i.e. regressions) are run
                     formula = lnoth.mort.p ~ age.for.pred + year) # To this formula
coefs2       <- ldply(regressions2, coef)                 # This captures the results of the 2 regressions in a dataset

#Apply to modelled dataset
modelled.mort$oth.mort2 <- 0                              # Place holder
modelled.mort$oth.mort2[modelled.mort$sex==1] <-          # Generating other mortality by sex
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==1],  # Make a 2424x3 matrix of 1s and ages and years
               modelled.mort$age[modelled.mort$sex==1],
               modelled.mort$year[modelled.mort$sex==1]), # 'exp' is to unlog model predictions
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[1,2],coefs2[1,3], coefs2[1,4]),  # Multiply by 3x1 matrix of sex
                ncol = 1, nrow = 3))                      # specific regression coefficients.
modelled.mort$oth.mort2[modelled.mort$sex==2] <-          
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==2],  
               modelled.mort$age[modelled.mort$sex==2],
               modelled.mort$year[modelled.mort$sex==2]), 
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[2,2],coefs2[2,3], coefs2[2,4]),  
                ncol = 1, nrow = 3))

# Check against year-specific model results
y <- 2012
par(mfrow=c(2,1))
plot(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==y],
     modelled.mort$oth.mort[modelled.mort$sex==1&modelled.mort$year==y],
     xlab = "Age", ylab = paste0("Other mortality rate, men, ",y), col = "green")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==y],
     modelled.mort$oth.mort2[modelled.mort$sex==1&modelled.mort$year==y],
     col = "red")
plot(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==y],
     modelled.mort$oth.mort[modelled.mort$sex==2&modelled.mort$year==y],
     xlab = "Age", ylab = paste0("Other mortality rate, women, ",y), col = "green")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==y],
       modelled.mort$oth.mort2[modelled.mort$sex==2&modelled.mort$year==y],
       col = "red") # And that is close enough.

# Check range of results if parameters vary by 99% confidence intervals. First, allowing alpha parameters
# only to vary.
# Get individual lm results by sex
lm.men   <- lm(lnoth.mort.p ~ age.for.pred + year, data = regress.db[regress.db$sex==1,])
lm.women <- lm(lnoth.mort.p ~ age.for.pred + year, data = regress.db[regress.db$sex==2,])

modelled.mort$oth.mort3 <- 0                              # Place holder for low confidence interval
modelled.mort$oth.mort3[modelled.mort$sex==1] <-          # Generating other mortality by sex
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==1],  # Make a 2424x3 matrix of 1s and ages and years
               modelled.mort$age[modelled.mort$sex==1],
               modelled.mort$year[modelled.mort$sex==1]), # 'exp' is to unlog model predictions
             ncol = 3, nrow = 2424)
      %*%matrix(c(as.numeric(lm.men$coefficients[1]-2.58*summary(lm.men)$coef[1,2]), 
                  coefs2[1,3], coefs2[1,4]),              # Multiply by 3x1 matrix of sex
                ncol = 1, nrow = 3))                      # specific regression coefficients.
modelled.mort$oth.mort3[modelled.mort$sex==2] <-          
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==2],  
               modelled.mort$age[modelled.mort$sex==2],
               modelled.mort$year[modelled.mort$sex==2]), 
             ncol = 3, nrow = 2424)
      %*%matrix(c(as.numeric(lm.women$coefficients[1]-2.58*summary(lm.women)$coef[1,2]), 
                  coefs2[2,3], coefs2[2,4]),  
                ncol = 1, nrow = 3))

modelled.mort$oth.mort4 <- 0                              # Place holder for high confidence interval
modelled.mort$oth.mort4[modelled.mort$sex==1] <-          # Generating other mortality by sex
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==1],  # Make a 2424x3 matrix of 1s and ages and years
               modelled.mort$age[modelled.mort$sex==1],
               modelled.mort$year[modelled.mort$sex==1]), # 'exp' is to unlog model predictions
             ncol = 3, nrow = 2424)
      %*%matrix(c(as.numeric(lm.men$coefficients[1]+2.58*summary(lm.men)$coef[1,2]),
                  coefs2[1,3], coefs2[1,4]),              # Multiply by 3x1 matrix of sex
                ncol = 1, nrow = 3))                      # specific regression coefficients.
modelled.mort$oth.mort4[modelled.mort$sex==2] <-          
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==2],  
               modelled.mort$age[modelled.mort$sex==2],
               modelled.mort$year[modelled.mort$sex==2]), 
             ncol = 3, nrow = 2424)
      %*%matrix(c(as.numeric(lm.women$coefficients[1]+2.58*summary(lm.women)$coef[1,2]),
                  coefs2[2,3], coefs2[2,4]),  
                ncol = 1, nrow = 3))

par(mfrow=c(2,1))
plot(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
     modelled.mort$oth.mort[modelled.mort$sex==1&modelled.mort$year==1998],
     xlab = "Age", ylab = "Other mortality rate, men, 1998", col = "green")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort2[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "red")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort3[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "blue")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort4[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "blue")
plot(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
     modelled.mort$oth.mort4[modelled.mort$sex==2&modelled.mort$year==1998],
     xlab = "Age", ylab = "Other mortality rate, women, 1998", col = "green")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort2[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "red")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort3[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "blue")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort4[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "blue")

# THIS SHOWS THAT THE RANGE IS TOTALLY LUDICROUS WHEN THE 99% COONFIDENCE INTERVALS ARE USED. INSTEAD GO 
# FOR +/- 30%

modelled.mort$oth.mort3 <- 0                              # Place holder for low confidence interval
modelled.mort$oth.mort3[modelled.mort$sex==1] <-          # Generating other mortality by sex
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==1],  # Make a 2424x3 matrix of 1s and ages and years
               modelled.mort$age[modelled.mort$sex==1],
               modelled.mort$year[modelled.mort$sex==1]), # 'exp' is to unlog model predictions
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[1,2]+log(0.7), coefs2[1,3], coefs2[1,4]), # Multiply by 3x1 matrix of sex
                ncol = 1, nrow = 3))                      # specific regression coefficients.
modelled.mort$oth.mort3[modelled.mort$sex==2] <-          
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==2],  
               modelled.mort$age[modelled.mort$sex==2],
               modelled.mort$year[modelled.mort$sex==2]), 
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[2,2]+log(0.7), coefs2[2,3], coefs2[2,4]),  
                ncol = 1, nrow = 3))

modelled.mort$oth.mort4 <- 0                              # Place holder for high confidence interval
modelled.mort$oth.mort4[modelled.mort$sex==1] <-          # Generating other mortality by sex
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==1],  # Make a 2424x3 matrix of 1s and ages and years
               modelled.mort$age[modelled.mort$sex==1],
               modelled.mort$year[modelled.mort$sex==1]), # 'exp' is to unlog model predictions
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[1,2]+log(1.3), coefs2[1,3], coefs2[1,4]),# Multiply by 3x1 matrix of sex
                ncol = 1, nrow = 3))                      # specific regression coefficients.
modelled.mort$oth.mort4[modelled.mort$sex==2] <-          
  exp(matrix(c(modelled.mort$cons[modelled.mort$sex==2],  
               modelled.mort$age[modelled.mort$sex==2],
               modelled.mort$year[modelled.mort$sex==2]), 
             ncol = 3, nrow = 2424)
      %*%matrix(c(coefs2[2,2]+log(1.3), coefs2[2,3], coefs2[2,4]),  
                ncol = 1, nrow = 3))

par(mfrow=c(2,1))
plot(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
     modelled.mort$oth.mort[modelled.mort$sex==1&modelled.mort$year==1998],
     xlab = "Age", ylab = "Other mortality rate, men, 1998", col = "green")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort2[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "red")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort3[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "blue")
points(modelled.mort$age[modelled.mort$sex==1&modelled.mort$year==1998],
       modelled.mort$oth.mort4[modelled.mort$sex==1&modelled.mort$year==1998],
       col = "blue")
plot(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
     modelled.mort$oth.mort[modelled.mort$sex==2&modelled.mort$year==1998],
     xlab = "Age", ylab = "Other mortality rate, women, 1998", col = "green")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort2[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "red")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort3[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "blue")
points(modelled.mort$age[modelled.mort$sex==2&modelled.mort$year==1998],
       modelled.mort$oth.mort4[modelled.mort$sex==2&modelled.mort$year==1998],
       col = "blue")


# Export parameters to a csv in Parameter space

alpha.mort.m <- c(coefs2[1,2]+log(0.7), coefs2[1,2]+log(1.3),  # Ln terms determine 30% either side of the best estimate.
                  NA, NA, NA, NA, NA, NA, NA, NA, "intercept")
alpha.mort.f <- c(coefs2[2,2]+log(0.7), coefs2[2,2]+log(1.3),  # Ln terms determine 30% either side of the best estimate.
                  NA, NA, NA, NA, NA, NA, NA, NA, "intercept")

# Combine into single dataset
parameter.space <- rbind(alpha.mort.m,
                         alpha.mort.f)

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

write.csv(parameter.space, "Other mortality parameters_10NOV2017.csv")


##################################################################################################
#                                                                                                #
# 2. Estimating migration rates                                                                  #
#                                                                                                #
##################################################################################################

# For this, we need to estimate single year estimates of total mortality by single year age and sex.
# This is done following the same method as was used for estimating the single year estimates of 'other mortality',
# so uses much of the same code, but based on 'mort' rather than 'oth.mort'

# Check plots of age against mortality rate (as percentage of population) in 1992 and 2015
age.r <- c(2.5,        # This is the age variable for regression analyses and plots. Under 1s ignored, since we already have single
           (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
           85)         # For the plots below, this is arbitrarily set as 85.

par(mfrow=c(2,2))
plot(age.r, 
     fert.mort$mort.p[(fert.mort$X%%5==3&fert.mort$year==1992&fert.mort$sex==1)|  # %% is modulo, takes every third reading
                      (fert.mort$age==85&fert.mort$year==1992&fert.mort$sex==1)], # This gets the final estimate for over 85s
     xlab = "Age", ylab = "Mort rate, 1992, men")
plot(age.r, 
     fert.mort$mort.p[(fert.mort$X%%5==3&fert.mort$year==1992&fert.mort$sex==2)|  
                      (fert.mort$age==85&fert.mort$year==1992&fert.mort$sex==2)],
     xlab = "Age", ylab = "Mort rate, 1992, women") 
plot(age.r, 
     fert.mort$mort.p[(fert.mort$X%%5==3&fert.mort$year==2015&fert.mort$sex==1)|  
                      (fert.mort$age==85&fert.mort$year==2015&fert.mort$sex==1)],
     xlab = "Age", ylab = "Mort rate, 2015, men") 
plot(age.r, 
     fert.mort$mort.p[(fert.mort$X%%5==3&fert.mort$year==2015&fert.mort$sex==2)|  
                      (fert.mort$age==85&fert.mort$year==2015&fert.mort$sex==2)], 
     xlab = "Age", ylab = "Mort rate, 2015, women") 

# Looks exponential in shape. Not clear what average age should be taken for the oldest group


# To find the parameter (y) for each year which minimises the square difference between predicted and measured rates:
sum.diff.men.by.j   <- c() # To store the results of the big loop initiated below
sum.diff.women.by.j <- c()

for (j in 65:95) {

y <- j                 # This sets the age for oldest age group

age.r <- c(2.5,        # This will be the age variable for regression analyses. Under 1s ignored, since we already have single
           (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
           y)          # See above.



# Run exponential fitting models
fert.mort$lnmort.p                        <- log(fert.mort$mort.p) # Produce log-transformed variable

fert.mort$age.for.pred                    <- fert.mort$age  # This introduces a new variable with age 85 replaced with y as 
fert.mort$age.for.pred[fert.mort$age==85] <- y              # average age for the over 85 group for the predicted rates.

fert.mort$pred.mort                       <- c(rep(0,4128)) # Dummy variable to be replaced with predicted rates in the loop below.
fert.mort$cons                            <- c(rep(1,4128)) # Column of 1s for the prediction process.


# Create reduced dataset for the regressions
regress.db  <- fert.mort[fert.mort$age%%5==3|fert.mort$age.for.pred==y,] # NB: last comma is to show I want all variables
regressions <- dlply(regress.db,                        # This runs a series of functions across diffent elements of regress.db
                     .(year, sex),                      # Split by year and sex
                     lm,                                # Linear models (i.e. regressions) are run
                     formula = lnmort.p ~ age.for.pred) # To this formula
coefs       <- ldply(regressions, coef)                 # This captures the results of the 48 regressions in a dataset

for (i in 1992:2015) {
  fert.mort$pred.mort[fert.mort$year==i&fert.mort$sex==1] <-                       # Replacing pred.mort by year and sex group
         exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==1],          # Make a 86x2 matrix of 1s and ages
                      fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==1]), # 'exp' is to unlog model predictions
                    ncol = 2, nrow = 86)
             %*%matrix(c(coefs[(2*(i-1992)+1),3],coefs[(2*(i-1992)+1),4]),         # Multiply by 2x1 matrix of year and sex
                       ncol = 1, nrow = 2))                                        # specific regression coefficients.
  
  fert.mort$pred.mort[fert.mort$year==i&fert.mort$sex==2] <-
    exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==2],
                 fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==2]),
               ncol = 2, nrow = 86)
        %*%matrix(c(coefs[(2*(i-1992)+2),3],coefs[(2*(i-1992)+2),4]),
                  ncol = 1, nrow = 2))
}

# Aim is to select y to minimise the difference between predicted and measured mortality rates up to 83, and then use that.
sum.diff.men        <- c() # To store the results of the small loop below
sum.diff.women      <- c()
for (i in 1992:2015){ # This is to calculate the sum of the square of differences between prediction and measure
  z1 <- sum(fert.mort$pred.mort[fert.mort$age%%5==3&fert.mort$sex==1&fert.mort$year==i]-
              fert.mort$mort.p[fert.mort$age%%5==3&fert.mort$sex==1&fert.mort$year==i])^2
  sum.diff.men <- append(sum.diff.men, z1)
  z2 <- sum(fert.mort$pred.mort[fert.mort$age%%5==3&fert.mort$sex==2&fert.mort$year==i]-
              fert.mort$mort.p[fert.mort$age%%5==3&fert.mort$sex==2&fert.mort$year==i])^2
  sum.diff.women <- append(sum.diff.women, z2)
}

sum.diff.men.by.j   <- rbind(sum.diff.men.by.j,   sum.diff.men)
sum.diff.women.by.j <- rbind(sum.diff.women.by.j, sum.diff.women)

}

sum.diff.men.by.j   <- as.data.frame(sum.diff.men.by.j)
sum.diff.women.by.j <- as.data.frame(sum.diff.women.by.j)

# Find minimum sum of squared difference for each year
min.men   <- c()
min.women <- c()
for (i in 1:24){
  z1 <- which.min(sum.diff.men.by.j[,i])+64
  min.men <- append(min.men, z1)
  z2 <- which.min(sum.diff.women.by.j[,i])+64
  min.women <- append(min.women, z2)
}

# Plot the results that show where the best exponential fit is
year <- 1992:2015
par(mfrow=c(1,2))
plot(year, min.men,   ylab = "Minimum squared diff by j, men")
plot(year, min.women, ylab = "Minimum squared diff by j, women")
# This demonstrates that for both men and women, the old age group should not be more than 85, but for many years you get a 
# better fit of the model if the oldest age group is less than 85. But that doesn't make sense, so from herein the oldest age is
# set as 85.

# Set the predicted mortality rates with oldest age group set as 85
y <- 85
age.r <- c(2.5,        # This will be the age variable for regression analyses. Under 1s ignored, since we already have single
           (1:16)*5+2, # year estimates. All other categories are set at mid-point of (mostly 5 year) age range.
           y)          # See above.
fert.mort$age.for.pred[fert.mort$age==85] <- y              # average age for the over 85 group for the predicted rates.
fert.mort$pred.mort                       <- c(rep(0,4128)) # Dummy variable to be replaced with predicted rates in the loop below.
fert.mort$cons                            <- c(rep(1,4128)) # Column of 1s for the prediction process.
# Create reduced dataset for the regressions
regress.db <- fert.mort[fert.mort$age%%5==3|fert.mort$age.for.pred==y,] # NB: last comma is to show I want all variables
regressions <- dlply(regress.db,                        # This runs a series of functions across diffent elements of regress.db
                     .(year, sex),                      # Split by year and sex
                     lm,                                # Linear models (i.e. regressions) are run
                     formula = lnmort.p ~ age.for.pred) # To this formula
coefs       <- ldply(regressions, coef)                 # This captures the results of the 48 regressions in a dataset

for (i in 1992:2015) {
  fert.mort$pred.mort[fert.mort$year==i&fert.mort$sex==1] <-                       # Replacing pred.mort by year and sex group
    exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==1],          # Make a 86x2 matrix of 1s and ages
                 fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==1]), # 'exp' is to unlog model predictions
               ncol = 2, nrow = 86)
        %*%matrix(c(coefs[(2*(i-1992)+1),3],coefs[(2*(i-1992)+1),4]),         # Multiply by 2x1 matrix of year and sex
                  ncol = 1, nrow = 2))                                        # specific regression coefficients.
  
  fert.mort$pred.mort[fert.mort$year==i&fert.mort$sex==2] <-
    exp(matrix(c(fert.mort$cons[fert.mort$year==i&fert.mort$sex==2],
                 fert.mort$age.for.pred[fert.mort$year==i&fert.mort$sex==2]),
               ncol = 2, nrow = 86)
        %*%matrix(c(coefs[(2*(i-1992)+2),3],coefs[(2*(i-1992)+2),4]),
                  ncol = 1, nrow = 2))
}

# Check predicted results against measured results for men and women in 1992 and 2015
# Red = measured results; Green = predicted results
par(mfrow=c(2,2))
plot(regress.db$age.for.pred[regress.db$year==1992&regress.db$sex==1], # Men, 1992
     regress.db$mort.p[regress.db$year==1992&regress.db$sex==1], 
     col = "red", xlab = "Age", ylab = "Mort rate, 1992, men")
points(fert.mort$age.for.pred[fert.mort$year==1992&fert.mort$sex==1],
       fert.mort$pred.mort[fert.mort$year==1992&fert.mort$sex==1],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==1992&regress.db$sex==2], # Women, 1992
     regress.db$mort.p[regress.db$year==1992&regress.db$sex==2], 
     col = "red", xlab = "Age", ylab = "Mort rate, 1992, women")
points(fert.mort$age.for.pred[fert.mort$year==1992&fert.mort$sex==2],
       fert.mort$pred.mort[fert.mort$year==1992&fert.mort$sex==2],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==2015&regress.db$sex==1], # Men, 2015
     regress.db$mort.p[regress.db$year==2015&regress.db$sex==1], 
     col = "red", xlab = "Age", ylab = "Mort rate, 2015, men")
points(fert.mort$age.for.pred[fert.mort$year==2015&fert.mort$sex==1],
       fert.mort$pred.mort[fert.mort$year==2015&fert.mort$sex==1],
       col = "green")
plot(regress.db$age.for.pred[regress.db$year==2015&regress.db$sex==2], # Women, 2015
     regress.db$mort.p[regress.db$year==2015&regress.db$sex==2], 
     col = "red", xlab = "Age", ylab = "Mort rate, 2015, women")
points(fert.mort$age.for.pred[fert.mort$year==2015&fert.mort$sex==2],
       fert.mort$pred.mort[fert.mort$year==2015&fert.mort$sex==2],
       col = "green")
# All close up to 83 but suggest that the oldest age mortality rate is not reached until much after 85 (especially in later years).
# For 2015 women, the rate reaches 0.15 at age 98.


## Estimate migration
# Let pop(x,y), mort(x,y) and mig(x,y) mean population, mortality and net migration for age group x and year y respectively.
# Then pop(x+1,y+1) = pop(x,y) + mig(x,y) - mort(x,y)
#          mig(x,y) = pop(x+1,y+1) - pop(x,y) + mort(x,y)

# First calculate absolute number of deaths by single year of age in dataset, from the modelled results above.
# NB: Use actual data for age 0 and 85 and over.
fert.mort$pred.mort.num <- 1 # Place holder
for (i in 1:4128) {
  if (fert.mort$age[i]==0 | fert.mort$age[i]==85) { fert.mort$pred.mort.num[i] <- fert.mort$mort[i] }
  else { fert.mort$pred.mort.num[i] <- round(fert.mort$population[i]*fert.mort$pred.mort[i], 0)}
}

# Generate migration estimates
fert.mort$mig <- 1 # Place holder
for (i in 1:4128) {
  if (fert.mort$year[i]==2015) { fert.mort$mig[i] <- NA }
  else if (fert.mort$age[i]==84 | fert.mort$age[i]==85) { fert.mort$mig[i] <- NA } # Can't use age 84 as 85 is not single year
  else { fert.mort$mig[i] <- fert.mort$population[86*2+1+i] - fert.mort$population[i] + fert.mort$pred.mort.num[i] }
}

# Plot migration by age, sex and year
par(mfrow=c(2,1))
plot(fert.mort$age[fert.mort$sex==1&fert.mort$year==1992], fert.mort$mig[fert.mort$sex==1&fert.mort$year==1992],
     xlab = "Age", ylab = "Migration, men", col=92, ylim = c(-20000,20000))
for (i in 1993:2014) {
  points(fert.mort$age[fert.mort$sex==1&fert.mort$year==i], fert.mort$mig[fert.mort$sex==1&fert.mort$year==i],
         col=i-1990)
}
plot(fert.mort$age[fert.mort$sex==2&fert.mort$year==1992], fert.mort$mig[fert.mort$sex==2&fert.mort$year==1992],
     xlab = "Age", ylab = "Migration, women", col=92, ylim = c(-20000,20000))
for (i in 1993:2014) {
  points(fert.mort$age[fert.mort$sex==2&fert.mort$year==i], fert.mort$mig[fert.mort$sex==2&fert.mort$year==i],
         col=i-1990)
}

# These migration figures can be used directly in the training period of the microsim.

##################################################################################################
#                                                                                                #
# 3. Estimating fertility rates                                                                  #
#                                                                                                #
##################################################################################################

# Really this is estimating how many new 18 year olds enter the model each year. These need to be
# estimated exactly for the training period and then  entered directly (i.e. not based on parameters
# that could vary) for the training period. Projecting forward doesn't matter as much.

# Produce plot of population of 18 year olds by year, 1992-2015
par(mfrow=c(2,1))
plot(1992:2015, fert.mort$population[fert.mort$X%%86==19&fert.mort$sex==1],
     xlab = "Year", ylab = "Population of 18 year olds, male")
plot(1992:2015, fert.mort$population[fert.mort$X%%86==19&fert.mort$sex==2],
     xlab = "Year", ylab = "Population of 18 year olds, female")