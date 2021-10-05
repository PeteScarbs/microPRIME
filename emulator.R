### emulator ###

# Last updated: 22nd April 2021

# This R script runs the emulator, which aims to identify a small set of vectors where the microsim
# produces results that are close to the results in the external dataset. The stages are described below.
#
# 0. If you want to play with previous emulated results, can just load them up here.
# 1. Generate a large sample of K vectors across the parameter space
# 2. Use the training data to produce scales for emulators of outcome variables
# 3. Run diagnostics of the outcome-specific emulators on the training dataset
# 4. Run the adequate emulators against the K vectors
# 5. Use the external dataset and estimates of uncertainty at different points in the process to estimate
#    the implausibility measure for each vector.
# 6. Produce graphs to explore how implausibility varies against parameters
# 7. Produce a subset of vectors where the implausibility measures for all outcomes with adequate
#    emulators are sufficiently small.
# 8. Produce a small sample from the subsample, which is then run through the microsim to generate new
#    training data.
# 9. Go back to stage 2 with the subsample of vectors and go around the loop until stopping criteria are
#    met.
#
# ADDITIONAL STAGE: If the subsample gets too small, then it may be necessary to generate a new sample
#                   by perturbing vectors that are included in the subsample.



##############################################################################
#                                                                            #
# Essential packages                                                         #
#                                                                            #
##############################################################################

# install.packages("emulator")
library("emulator")
# install.packages("lhs")
library("lhs")
# install.packages("plyr")
library("plyr")
# install.packages("dplyr")
library("dplyr")


##############################################################################
#                                                                            #
# Initial settings                                                           #
#                                                                            #
##############################################################################

work <- "desktop"# Choice between 'desktop' and 'laptop'. Ensures script uses correct directories.
new  <- "no"    # If "yes" then this starts a whole new loop and generates a new sample of K points to emulate at
                 # If "no" then this uses saved points from a previous loop
K    <- 1000000  # This is the size of the initial sample, which will be reduced in each of the waves
                 # of the cycle.
p    <- 19       # Number of varying parameters included in the microsim per gender (i.e. number of variables in train.data.input.m/f)
n.m  <- 22       # The number of outputs estimated by the training data for men 
n.f  <- 22       # The number of outputs estimated by the training data for women
n.train <- 450   # The number of runs of the model in the training data
r.min <- 0.6     # The minimum acceptable correlation between the emulated data and the training data
imp.max <- 3     # The maximum value of implausibility alowed for any single output. NB, Andrianakis et
                 # al recommend this value should be 3.
wave.number <- 2 # Set the wave number for saving the diagnostic images
model.number<- 18# Sets the model number (see Model Run Notes) for saving images
parameter.space.m <- "parameters_m_12MAR2021_worst.csv" # The name of the csv file produced by parameter draw.R
                                                  # that contains the parameter space for men
parameter.space.f <- "parameters_f_12MAR2021_worst.csv" # The name of the csv file produced by parameter draw.R
                                                  # that contains the parameter space for women
train.data.input.m  <- "train.input.m_17APR2021.csv" # The name of the csv which contains the specific parameters
                                                   # used in the training data for men
train.data.input.f  <- "train.input.f_17APR2021.csv" # The name of the csv which contains the specific parameters
                                                   # used in the training data for women
# See training data set up.R for how these following four files are constructed from microsim output.
train.data.output.m <- "train.output.m_17APR2021.csv" # The name of the csv that contains the
                                                 # outcomes of the training data for men
train.data.output.f <- "train.output.f_17APR2021.csv" # The name of the csv that contains the
                                                 # outcomes of the training data for women
train.data.var.m <- "train.output.var.m_17APR2021.csv" # The name of the csv that contains the
                                                 # variance of the outcomes of the training data for men
train.data.var.f <- "train.output.var.f_17APR2021.csv" # The name of the csv that contains the
                                                 # variance of the outcomes of the training data for women
external.data.m <- "external men_12DEC2019.csv"  # The name of the csv that contains the external dataset
                                                 # estimates of outputs for men, including estimates of error
external.data.f <- "external women_12DEC2019.csv" # The name of the csv that contains the external dataset
                                                 # estimates of outputs for women, including estimates of error
use.saved.emulators <- "No"                      # Set to "Yes" if want to skip the emulation section, "No" otherwise.
saved.em.m      <- "Emulator men_30MAR2020.csv"   # File name of saved emulators to be used, men.
saved.em.f      <- "Emulator women_30MAR2020.csv" # File name of saved emulators to be used, women.
saved.r.m       <- "Correlation data men_30MAR2020.csv" # File name for saved diagnostics.
saved.r.f       <- "Correlation data women_30MAR2020.csv"
saved.scales.m  <- "Scales men_30MAR2020.csv"     # File name for saved scales
saved.scales.f  <- "Scales women_30MAR2020.csv"
imp.m           <- "Implausible men_24JUL2019.csv" # File name for saved implausibility measures
imp.f           <- "Implausible women_24JUL2019.csv"
remaining.vectors.m <- "Remaining vector space_m_09APR2021.csv" # These are the name of the output from the previous
remaining.vectors.f <- "Remaining vector space_f_09APR2021.csv" # wave of emulator, which sets the remaining vector
                                                                # for this wave.
# Note, to set the names of the output (remaining vector space and random draws) change lines 422, 523 and 555)


##############################################################################
#                                                                            #
# 0. Load previous emulated results, diagnostics and implausibility data     #
#                                                                            #
##############################################################################

if (use.saved.emulators=="Yes") {
  if (work=="desktop") { # This sets the directory according to the initial settings
    setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  } else {
    if (work=="laptop"){
      setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
    }  else { "Error" }
  }
  Em.f <- read.csv(saved.em.f)
  Em.f <- Em.f[,2:(n.f+1)] # This removes the first column, which is not needed.
  Em.m <- read.csv(saved.em.m)
  Em.m <- Em.m[,2:(n.m+1)] # This removes the first column, which is not needed.
  implausible.m <- read.csv(imp.m)
  implausible.m <- implausible.m[,2:(n.m+2)] # Removes first column. Needs '+2' as there is the max variable.
  implausible.f <- read.csv(imp.f)
  implausible.f <- implausible.f[,2:(n.f+2)]
  diagnostics.m <- read.csv(saved.r.m)
  diagnostics.f <- read.csv(saved.r.f)
  scales.m      <- read.csv(saved.scales.m)
  scales.f      <- read.csv(saved.scales.f)
  } else if (use.saved.emulators=="No") {
    print("This section is skipped")
  } else { 
    print("Error in use.saved.emulators variable") 
    }
 

##############################################################################
#                                                                            #
# 1. Generate a sample of K vectors across the parameter space               #
#                                                                            #
##############################################################################

if (new=="yes") { # This states that the code below will only run if new is set to "yes"

X.m <- randomLHS(K,p)
X.f <- randomLHS(K,p) # Generates two latin hypercubes for men and women with K points over 
                      # p dimensions of [0,1] space
## Apply X to the parameter space.
# Load the parameter space dataset
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Parameter space")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Parameter space")
  }  else { "Error" }
}
parameters.m           <- read.csv(parameter.space.m)
parameters.f           <- read.csv(parameter.space.f)
names                  <- c("N", "name", "low", "high")
colnames(parameters.m) <- names
colnames(parameters.f) <- names

# Translate X from [0,1] to parameter space
for (j in 1:p) {
  X.m[,j] <- parameters.m[j,3] + X.m[,j]*(parameters.m[j,4]-parameters.m[j,3])
  X.f[,j] <- parameters.f[j,3] + X.f[,j]*(parameters.f[j,4]-parameters.f[j,3])
}

} else if (new=="no") {
  if (work=="desktop") { # This sets the directory according to the initial settings
    setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Training data")
  } else {
    if (work=="laptop"){
      setwd("F:\\microPRIME\\Data repositories\\Training data")
    }  else { "Error" }
  }
  
  v.remain.m <- read.csv(remaining.vectors.m)
  v.remain.f <- read.csv(remaining.vectors.f)
  v.remain.m <- v.remain.m[,2:(p+1)]
  v.remain.f <- v.remain.f[,2:(p+1)]
  
} else {
  print("Error in variable new")
}

##############################################################################
#                                                                            #
# 2. Produce scales for emulators of outcome variables                       #
#                                                                            #
##############################################################################

# Set up the inputs and outputs from the training data
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Training data")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Training data")
  }  else { "Error" }
}

train.input.m  <- read.csv(train.data.input.m)
train.input.f  <- read.csv(train.data.input.f)
train.output.m <- read.csv(train.data.output.m)
train.output.f <- read.csv(train.data.output.f) 

# Remove redundant first variable from each
train.input.m  <- train.input.m[,2:(p+1)]
train.input.f  <- train.input.f[,2:(p+1)]
train.output.m <- train.output.m[,2:(n.m+1)]
train.output.f <- train.output.f[,2:(n.f+1)]

# Limit the input data to the number of iterations in the most recent run
train.input.m  <- train.input.m[1:n.train,]
train.input.f  <- train.input.f[1:n.train,]

# Generate scales for each outcome

##############################################################################
# Note that the number of iterations selected for the scale finder can influence
# how well the emulated results fit the observed results, with greater number
# of iterations generally resulting in a better emulated fit, but taking longer.
# Set the number of iterations here

# Set timer for emulator run
model.run.start <- Sys.time()

# Scales for men
scales.it.m <- 100

# Scales for women
scales.it.f <- 100

# Set up dataset to collect results
scales.m <- NULL
scales.f <- NULL

for (i in 1:n.m) {
  parscale <- optimal.scales(val=as.matrix(train.input.m), scales.start=rep(1,p), d=train.output.m[,i],
                                       method="SANN",control=list(trace=1000,maxit=scales.it.m), give=FALSE)
  scales.m <- cbind(scales.m, parscale)
  scales.m <- as.data.frame(scales.m)
  }
for (i in 1:n.f) {
  parscale <- optimal.scales(val=as.matrix(train.input.f), scales.start=rep(1,p), d=train.output.f[,i], 
                                       method="SANN",control=list(trace=1000,maxit=scales.it.f), give=FALSE)
  scales.f <- cbind(scales.f, parscale)
  scales.f <- as.data.frame(scales.f)
}

# Save output
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  }  else { "Error" }
}
write.csv(scales.m, "Scales men_17APR2021.csv")
write.csv(scales.f, "Scales women_17APR2021.csv")

##############################################################################
#                                                                            #
# 2. Diagnostics of emulators against training data                          #
#                                                                            #
##############################################################################

# Assess agreement between emulated points and actual evaluations from training data.
# Code below returns Pearson's r, Spearman's rho and Kendall's tau for each output variable.
# Note that Spearman and Kendall are non-parametric and so not influenced by outliers.
# At present, remaining code is based on Pearson, so Kendall and Spearmen provided purely for observation.

# Set up dataframes to collect results
diagnostics.m <- data.frame(outcome = 1:n.m,
                            r = numeric(n.m),
                            rho = numeric(n.m),
                            tau = numeric(n.m))

diagnostics.f <- data.frame(outcome = 1:n.f,
                            r = numeric(n.f),
                            rho = numeric(n.f),
                            tau = numeric(n.f))

for (i in 1:n.m) {
  assign(paste0("A.m",i), corr.matrix(as.matrix(train.input.m), scales=scales.m[,i]))
  assign(paste0("d.m",i), estimator(as.matrix(train.input.m), get(paste0("A.m",i)), train.output.m[,i],
                                    scales=scales.m[,i]))
  r <- cor(cbind(get(paste0("d.m",i)),train.output.m[,i]), method = "pearson")
  diagnostics.m$r[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Pearson r for variable ",i,":",round(r[1,2],3)))
  r <- cor(cbind(get(paste0("d.m",i)),train.output.m[,i]), method = "spearman")
  diagnostics.m$rho[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Spearman rho for variable ",i,":",round(r[1,2],3)))
  r <- cor(cbind(get(paste0("d.m",i)),train.output.m[,i]), method = "kendall")
  diagnostics.m$tau[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Kendall tau for variable ",i,":",round(r[1,2],3)))
}

for (i in 1:n.f) {
  assign(paste0("A.f",i), corr.matrix(as.matrix(train.input.f), scales=scales.f[,i]))
  assign(paste0("d.f",i), estimator(as.matrix(train.input.f), get(paste0("A.f",i)), train.output.f[,i], 
                                    scales=scales.f[,i]))
  r <- cor(cbind(get(paste0("d.f",i)),train.output.f[,i]), method = "pearson")
  diagnostics.f$r[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Pearson r for variable ",i,":",round(r[1,2],3)))
  r <- cor(cbind(get(paste0("d.f",i)),train.output.f[,i]), method = "spearman")
  diagnostics.f$rho[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Spearman rho for variable ",i,":",round(r[1,2],3)))
  r <- cor(cbind(get(paste0("d.f",i)),train.output.f[,i]), method = "kendall")
  diagnostics.f$tau[i] <- r[1,2] # This saves the r value for each outcome.
  print(paste0("Kendall tau for variable ",i,":",round(r[1,2],3)))
}

# Save output
write.csv(diagnostics.m, file = "Correlation data men_17APR2021.csv")
write.csv(diagnostics.f, file = "Correlation data women_17APR2021.csv")

# Produce images of the fit for the model run notes and save them
for(i in 1:n.m){
  if (diagnostics.m$r[i]<r.min) { # i.e. if fit is inadequate
    print(paste0("Inadequate fit for outcome ",i))
   }
  else {
    if (work=="desktop"){
      mypath <- file.path("K:","CPNP","Pete","MODELLING ALL RISK FACTORS","microPRIME","Model run notes","IMAGES", paste0("MODEL ",model.number), # Sets the file path
                          paste0("Wave ", wave.number, ", Variable ", i, ", Men.jpg")) # Sets the name of the saved image
    } else if (work=="laptop") {
      mypath <- file.path("F:","microPRIME","Model run notes","IMAGES", paste0("MODEL ",model.number), # Sets the file path
                          paste0("Wave ", wave.number, ", Variable ", i, ", Men.jpg")) # Sets the name of the saved image
    } else {
      print("Error in filename")
    }
    jpeg(file=mypath)
    mytitle = paste0("Wave ", wave.number, ", Variable ", i, ", Men, r = ",round(diagnostics.m$r[i],2))
    plot(get(paste0("d.m",i)), 
         train.output.m[,i], 
         xlab = "Emulated results", ylab = "Microsimulated results", 
         main = mytitle)
    dev.off()
  }
}

for(i in 1:n.f){
  if (diagnostics.f$r[i]<r.min) { # i.e. if fit is inadequate
    print(paste0("Inadequate fit for outcome ",i))
  }
  else {
    if (work=="desktop"){
      mypath <- file.path("K:","CPNP","Pete","MODELLING ALL RISK FACTORS","microPRIME","Model run notes","IMAGES", paste0("MODEL ",model.number), # Sets the file path
                          paste0("Wave ", wave.number, ", Variable ", i, ", Women.jpg")) # Sets the name of the saved image
    } else if (work=="laptop") {
      mypath <- file.path("F:","microPRIME","Model run notes","IMAGES", paste0("MODEL ",model.number), # Sets the file path
                          paste0("Wave ", wave.number, ", Variable ", i, ", Women.jpg")) # Sets the name of the saved image
    } else {
      print("Error in filename")
    }
    jpeg(file=mypath)
    mytitle = paste0("Wave ", wave.number, ", Variable ", i, ", Women, r = ",round(diagnostics.f$r[i],2))
    plot(get(paste0("d.f",i)), 
         train.output.f[,i], 
         xlab = "Emulated results", ylab = "Microsimulated results", 
         main = mytitle)
    dev.off()
  }
}

##############################################################################
#                                                                            #
# 4. Run emulators for adequate outcomes against the remaining vectors       #
#                                                                            #
##############################################################################

Sys.time()
if (use.saved.emulators == "Yes") { # If Yes this will load saved results rather than
                                    # run the emulators to save time.
  if (work=="desktop") { # This sets the directory according to the initial settings
    setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  } else {
    if (work=="laptop"){
      setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
    }  else { "Error" }
  }
  Em.m <- read.csv(saved.em.m)
  Em.m <- Em.m[,2:(n.m+1)] # This removes the first column, which is not needed.
} else if (new=="yes") {
  Em.m <- NULL  # Set up a folder to save to.
for (i in 1:n.m) {
  if (diagnostics.m$r[i] < r.min) { # i.e. if fit is inadequate
    print(paste0("Inadequate fit for outcome ",i))
    x <- c(rep(99999,K)) # Fill a vector with 99999 (doesn't work with NULL)
    }
  else {
    x <- interpolant.quick(X.m,                                 # The set of points that the emulator is to be run at
                           train.output.m[,i],                  # The vector of observations
                           train.input.m,                       # The full set of measured inputs
                           Ainv = solve(get(paste0("A.m",i))),  # The covariance matrix
                           scales = scales.m[,i]) # The scales
    print(i)
  }
  Em.m <- cbind(Em.m,x) # Combine results
}

} else if(new=="no") {
  Em.m <- NULL  # Set up a folder to save to.
  for (i in 1:n.m) {
    if (diagnostics.m$r[i] < r.min) { # i.e. if fit is inadequate
      print(paste0("Inadequate fit for outcome ",i))
      x <- c(rep(99999,length(v.remain.m[,1]))) # Fill a vector with 99999 (doesn't work with NULL)
    }
    else {
      x <- interpolant.quick(as.matrix(v.remain.m),               # The set of points that the emulator is to be run at
                             train.output.m[,i],                  # The vector of observations
                             train.input.m,                       # The full set of measured inputs
                             Ainv = solve(get(paste0("A.m",i))),  # The covariance matrix
                             scales = scales.m[,i]) # The scales
      print(i)
    }
    Em.m <- cbind(Em.m,x) # Combine results
  }
} else {
  print("Error in variable new")
}
  
if (use.saved.emulators == "Yes") { # If Yes this will load saved results rather than
  # run the emulators to save time.
  if (work=="desktop") { # This sets the directory according to the initial settings
    setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  } else {
    if (work=="laptop"){
      setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
    }  else { "Error" }
  }
  Em.f <- read.csv(saved.em.f)
  Em.f <- Em.f[,2:(n.f+1)] # This removes the first column, which is not needed.
} else if (new=="yes") {
  Em.f <- NULL  # Set up a folder to save to.
  for (i in 1:n.f) {
    if (diagnostics.f$r[i] < r.min) { # i.e. if fit is inadequate
      print(paste0("Inadequate fit for outcome ",i))
      x <- c(rep(99999,K)) # Fill a vector with 99999 (doesn't work with NULL)
    }
    else {
      x <- interpolant.quick(X.f,                                 # The set of points that the emulator is to be run at
                             train.output.f[,i],                  # The vector of observations
                             train.input.f,                       # The full set of measured inputs
                             Ainv = solve(get(paste0("A.f",i))),  # The covariance matrix
                             scales = scales.f[,i]) # The scales
      print(i)
    }
    Em.f <- cbind(Em.f,x) # Combine results
  }
  
} else if (new=="no") {
  Em.f <- NULL  # Set up a folder to save to.
  for (i in 1:n.f) {
    if (diagnostics.f$r[i] < r.min) { # i.e. if fit is inadequate
      print(paste0("Inadequate fit for outcome ",i))
      x <- c(rep(99999,length(v.remain.f[,1]))) # Fill a vector with 99999 (doesn't work with NULL)
    }
    else {
      x <- interpolant.quick(as.matrix(v.remain.f),               # The set of points that the emulator is to be run at
                             train.output.f[,i],                  # The vector of observations
                             train.input.f,                       # The full set of measured inputs
                             Ainv = solve(get(paste0("A.f",i))),  # The covariance matrix
                             scales = scales.f[,i]) # The scales
      print(i)
    }
    Em.f <- cbind(Em.f,x) # Combine results
  }
} else {
  print("Error in variable new")
}

Sys.time() # For one million vectors, this takes about half an hour per outcome.

# Save output
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  }  else { "Error" }
}

write.csv(Em.m, file = "Emulator men_17APR2021.csv")
write.csv(Em.f, file = "Emulator women_17APR2021.csv")

##############################################################################
#                                                                            #
# 5. Estimate implausibility measures                                        #
#                                                                            #
##############################################################################

## Process:
# Need to identify the standard deviation for each of the following:
# 1. Observation uncertainty (i.e. uncertainty in the external dataset)
# 2. Code uncertainty (i.e. uncertainty in the emulated variables)
# 3. Ensemble uncertainty (i.e. uncertainty in the training data)
# Once those have been estimated, they are combined in the implausibility measure and then
# estimates are deemed implausible if the measure is greater than imp.max (usually set as 
# 3, since there is a 95% probability that a real result will be less than 3 standard 
# deviations away - see Andrianakis)

# Estimate observation uncertainty
# Open dataframe of external data
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\External dataset")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\External dataset")
  }  else { "Error" }
}
external.m <- read.csv(external.data.m)
external.f <- read.csv(external.data.f) # Note that for outcome i, the observation uncertainty
                                        # is external.m[i,4]

# Estimate code uncertainty for all outputs (i.e. uncertainty in the emulator)
# As shown in the emulator examples.R script (lines 1024-1029), the variance is the same at all points
# that are emulated, so it is just necessary to estimate the variance at a single point.

codeu.m <- NULL
for (i in 1:n.m) {
  if (diagnostics.m$r[i] < r.min) {
    codeu.m <- rbind(codeu.m, 99999) # Note that a number is used so that the variance estimates are
                                     # stored as numbers.
  }
  else {
    a <- interpolant(as.vector(unlist(train.input.m[1,])),  # The vector of parameters where I want an emulated result
                     train.output.m[,i],    # The vector of observations 
                     matrix(unlist(train.input.m), nrow = n.train, ncol = p), # The full set of training data
                     Ainv = solve(get(paste0("A.m",i))),  # The covariance matrix
                     scales = scales.m[,i], # The scales, 
                     give.full.list = TRUE) # Needed to get the variance estimate
    codeu.m <- rbind(codeu.m,a$sigmahat.square)
  }
}

codeu.f <- NULL
for (i in 1:n.f) {
  if (diagnostics.f$r[i] < r.min) {
    codeu.f <- rbind(codeu.f, 99999) # Note that a number is used so that the variance estimates are
    # stored as numbers.
  }
  else {
    a <- interpolant(as.vector(unlist(train.input.f[1,])),     # The vector of parameters where I want an emulated result
                     train.output.f[,i],    # The vector of observations 
                     matrix(unlist(train.input.f), nrow = n.train, ncol = p), # The full set of training data
                     Ainv = solve(get(paste0("A.f",i))),  # The covariance matrix
                     scales = scales.f[,i],
                     give.full.list = TRUE) # Needed to get the variance estimate
    codeu.f <- rbind(codeu.f,a$sigmahat.square)
  }
} 

## Estimate ensemble uncertainty (i.e. uncertainty in the training data)
# In Andrianakis et al, for each vector of training inputs, they run the model many times to account for
# stochastic variability. The ensemble uncertainty is estimated to be the variance in the outputs from
# the many model runs. However, my model doesn't have any interactions between agents, so in effect
# running one model with 50,000 agents is identical to running 5 models with 10,000 agents each. For that
# reason, it makes more sense for me to measure this as uncertainty in the training data outputs from
# a single run, using standard methods to account for sample variability.

# Open datasets

if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Training data")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Training data")
  }  else { "Error" }
}

var.m <- read.csv(train.data.var.m)
var.f <- read.csv(train.data.var.f) # The ensemble uncertainty can be taken as the mean standard deviation
                                    # across all of the training data
var.m <- var.m[,2:(n.m+1)]
var.f <- var.f[,2:(n.f+1)] # Get rid of useless first variable


# Calculate implausibility measure across all emulated datapoints
# I(x) = abs(z - g(x))/((var(c(x))+var(e(x))+var(o(x)))^0.5), where:
# x is a vector of all parameters used in the exercise
# z is the estimate from the external dataset that we are trying to replicate
# c(x) is the code uncertainty
# e(x) is the ensemble uncertainty
# o(x) is the observation uncertainty

implausible.m <- NULL
implausible.f <- NULL

for (i in 1:n.m){
  if (diagnostics.m$r[i] < r.min) {
    if (new=="yes") {
      y <- rep(-99999,K) # NB: These are set as -99999 so that they are never the maximum figure when combining
                       # implausibility estimates across outcomes. See below
      implausible.m <- cbind(implausible.m,y)
    } else if (new=="no") {
      y <- rep(-99999,length(v.remain.m[,1])) # NB: These are set as -99999 so that they are never the maximum figure when combining
      # implausibility estimates across outcomes. See below
      implausible.m <- cbind(implausible.m,y)
    } else {
      print("Error in variable new")
    }
  }
  else {
    y <- abs(external.m[i,3]-Em.m[,i])/((codeu.m[i]+(mean(var.m[,i]))+external.m[i,4])^0.5)
    implausible.m <- cbind(implausible.m,y)
  }
}

for (i in 1:n.f){
  if (diagnostics.f$r[i] < r.min) {
    if (new=="yes") {
      y <- rep(-99999,K)
      implausible.f <- cbind(implausible.f,y)
    } else if (new=="no") {
      y <- rep(-99999,length(v.remain.f[,1]))
      implausible.f <- cbind(implausible.f,y)
    } else {
      print("Error in variable new")
    }
  }
  else {
    y <- abs(external.f[i,3]-Em.f[,i])/((codeu.f[i]+(mean(var.f[,i]))+external.f[i,4])^0.5)
    implausible.f <- cbind(implausible.f,y)
  }
}

names.m <- c(paste0("v",1:n.m))
names.f <- c(paste0("v",1:n.f))

implausible.m <- data.frame(implausible.m)
implausible.f <- data.frame(implausible.f)

colnames(implausible.m) <- names.m
colnames(implausible.f) <- names.f

# Determine subsample of emulated points where implausibility is sufficiently small for all outputs
# First find maximum implausibility across all outcomes.

implausible.m$max <- apply(implausible.m, # The dataset over which to apply the function
                           1,             # Indicates that function should be aplpied over rows (2 would be columns)
                           max)           # The function to be applied
implausible.f$max <- apply(implausible.f,1,max)

# Save output
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\R modules\\Emulator\\Saved emulator results")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\R modules\\Emulator\\Saved emulator results")
  }  else { "Error" }
}
write.csv(implausible.m, file = "Implausible men_17APR2021.csv")
write.csv(implausible.f, file = "Implausible women_17APR2021.csv")

# End timer
model.run.end <- Sys.time() 
model.run.end - model.run.start

##############################################################################
#                                                                            #
# 7. Produce subset of vectors where implausibility is sufficiently small    #
#                                                                            #
##############################################################################

# NB: THIS CODE BELOW DOES NOT WORK
if (new=="yes"){
  ifelse(as.numeric(quantile(implausible.m$max, prob = 0.1))<imp.max, # If 10% threshold is less than imp.max threshold
         X.m.sub <- X.m[implausible.m$max<imp.max,], # Then take imp.max threshold
         X.m.sub <- X.m[implausible.m$max<=as.numeric(quantile(implausible.m$max, prob = 0.1)),]) # o/w take 10% threshold
  ifelse(as.numeric(quantile(implausible.f$max, prob = 0.1))<imp.max,
         X.f.sub <- X.f[implausible.f$max<imp.max,],
         X.f.sub <- X.f[implausible.f$max<=as.numeric(quantile(implausible.f$max, prob = 0.1)),])
} else if (new=="no") {
  ifelse(as.numeric(quantile(implausible.m$max, prob = 0.1))<imp.max,
         X.m.sub <- v.remain.m[implausible.m$max<imp.max,],
         X.m.sub <- v.remain.m[implausible.m$max<=as.numeric(quantile(implausible.m$max, prob = 0.1)),])
  ifelse(as.numeric(quantile(implausible.f$max, prob = 0.1))<imp.max,
         X.f.sub <- v.remain.f[implausible.f$max<imp.max,],
         X.f.sub <- v.remain.f[implausible.f$max<=as.numeric(quantile(implausible.f$max, prob = 0.1)),])
} else {
  print("Problem with variable new")
}

# Save these for following iteration
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Training data")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Training data")
  }  else { "Error" }
}

write.csv(X.m.sub, file = "Remaining vector space_m_17APR2021.csv")
write.csv(X.f.sub, file = "Remaining vector space_f_17APR2021.csv")

##############################################################################
#                                                                            #
# 8. Produce subsample of vectors for next training data                     #
#                                                                            #
##############################################################################

sample.m <- sample_n(as.data.frame(X.m.sub), # The data frame to be sampled from
                     450,                    # The number of sampled rows - should be at least 10 times the number of parameters
                     replace = FALSE)        # This means that a row cannot be sampled twice

sample.f <- sample_n(as.data.frame(X.f.sub),450,replace = FALSE)

# Name variables
parameter.names.m <- c("alpha.30dcf.m",
                       #"alpha.mort.m",
                       #"treatment.inc.mi.m",
                       "treatment.cf.mi.m",
                       "bmi.m.r1525",
                       "bmi.m.r2550",
                       "sbp.alpha.m",
                       "tcho.alpha.m",
                       "smok.alpha.m",
                       "diab.m",
                       "bmi.trend.age.m",
                       "bmi.trend.year.m",
                       #"bmi.sdtrend.age.m",
                       #"bmi.sdtrend.year.m",
                       "sbp.trend.age.m",
                       "sbp.trend.year.m",
                       "tcho.trend.age.m",
                       "tcho.trend.year.m",
                       "smok.trend.age.m",
                       "smok.trend.year.m",
                       "diab.trend.age.m",
                       "diab.trend.year.m",
                       "prev.alpha.m")
parameter.names.f <- c("alpha.30dcf.f",
                       #"alpha.mort.f",
                       #"treatment.inc.mi.f",
                       "treatment.cf.mi.f",
                       "bmi.f.r1525",
                       "bmi.f.r2550",
                       "sbp.alpha.f",
                       "tcho.alpha.f",
                       "smok.alpha.f",
                       "diab.f",
                       "bmi.trend.age.f",
                       "bmi.trend.year.f",
                       #"bmi.sdtrend.age.f",
                       #"bmi.sdtrend.year.f",
                       "sbp.trend.age.f",
                       "sbp.trend.year.f",
                       "tcho.trend.age.f",
                       "tcho.trend.year.f",
                       "smok.trend.age.f",
                       "smok.trend.year.f",
                       "diab.trend.age.f",
                       "diab.trend.year.f",
                       "prev.alpha.f")
colnames(sample.m) <- parameter.names.m
colnames(sample.f) <- parameter.names.f

# Combine
sample <- cbind(sample.m, sample.f)

# Export to csv
if (work=="desktop") { # This sets the directory according to the initial settings
  setwd("K:\\CPNP\\Pete\\MODELLING ALL RISK FACTORS\\microPRIME\\Data repositories\\Random Draws")
} else {
  if (work=="laptop"){
    setwd("F:\\microPRIME\\Data repositories\\Random Draws")
  }  else { "Error" }
}

write.csv(sample, file = "Random draws_22APR2021.csv")

length(X.m.sub[,1])
length(X.f.sub[,1]) # To collect the number of remaining vectors

quantile(implausible.m$max, prob = 0.1) # To collect the threshold for 10% of remaining vectors
quantile(implausible.f$max, prob = 0.1)




##############################################################################
#                                                                            #
# Produce graphs of emulated outcomes versus external datasets               #
#                                                                            #
##############################################################################

title.name <- c("1993, 55-64",
                "1993, 65-74",
                "1993, 75-84",
                "1993, 85+",
                "1994, 55-64",
                "1994, 65-74",
                "1994, 75-84",
                "1994, 85+",
                "1998, 55-64",
                "1998, 65-74",
                "1998, 75-84",
                "1998, 85+",
                "2003, 55-64",
                "2003, 65-74",
                "2003, 75-84",
                "2003, 85+",
                "2005, 65-74",
                "2005, 75-84",
                "2005, 85+",
                "2006, 55-64",
                "2006, 65-74",
                "2006, 75-84",
                "2006, 85+",
                "2011, 55-64",
                "2011, 65-74",
                "2011, 75-84",
                "2011, 85+")

par(mfrow = c(7,4),
    mar = c(1,1,1,1))

for (i in 1:16){
  hist(Em.m[,i+4], xlim = c(0,1), xaxt = "n", main = title.name[i], xlab = NULL, ylab = NULL)
  axis(1, at=0:1)
  abline(v = external.m[i+4,3], col = "Red")
}

hist(implausible.m$max, main = "Maximum implausibility", xaxt = "n")

for (i in 17:27) {
  hist(Em.m[,i+4], xlim = c(0,1), xaxt = "n", main = title.name[i], xlab = NULL, ylab = NULL)
  axis(1, at=0:1)
  abline(v = external.m[i+4,3], col = "Red")
}

par(mfrow = c(7,4),
    mar = c(1,1,1,1))

for (i in 1:16){
  hist(Em.f[,i+4], xlim = c(0,1), xaxt = "n", main = title.name[i], xlab = NULL, ylab = NULL)
  axis(1, at=0:1)
  abline(v = external.f[i+4,3], col = "Red")
}

hist(implausible.f$max, main = "Maximum implausibility", xaxt = "n")

for (i in 17:27) {
  hist(Em.f[,i+4], xlim = c(0,1), xaxt = "n", main = title.name[i], xlab = NULL, ylab = NULL)
  axis(1, at=0:1)
  abline(v = external.f[i+4,3], col = "Red")
}


##############################################################################
#                                                                            #
# Produce histograms of implausibility and emulated outcomes                 #
#                                                                            #
##############################################################################

for (i in 1:n.m) {
  if (diagnostics.m$r[i]<r.min) {
    print(paste0("No emulated results for outcome",i))
  }
  else if (external.m[i,5] == "Prevalence") {
    hist(Em.m[,i], 
         main = paste0(external.m[i,5], ", ", external.m[i,1], ", ", external.m[i,2]),
         xlim = c(0,0.4), col = "Light blue")
    abline(v = external.m[i,3], col = "Red")
    abline(v = (external.m[i,3]-1.96*external.m[i,4]^0.5), col = "Red", lty = 2)
    abline(v = (external.m[i,3]+1.96*external.m[i,4]^0.5), col = "Red", lty = 2)
  }
  else {
    hist(Em.m[,i], 
         main = paste0(external.m[i,5], ", ", external.m[i,1], ", ", external.m[i,2]),
         col = "Light blue")
    abline(v = external.m[i,3], col = "Red")
    abline(v = (external.m[i,3]-1.96*external.m[i,4]^0.5), col = "Red", lty = 2)
    abline(v = (external.m[i,3]+1.96*external.m[i,4]^0.5), col = "Red", lty = 2)
    
  }
}

for (i in 1:n.f) {
  if (diagnostics.f$r[i]<r.min) {
    print(paste0("No emulated results for outcome",i))
  }
  else if (external.f[i,5] == "Prevalence") {
    hist(Em.f[,i], 
         main = paste0(external.f[i,5], ", ", external.f[i,1], ", ", external.f[i,2]),
         xlim = c(0,0.4), col = "Pink")
    abline(v = external.f[i,3], col = "Red")
    abline(v = (external.f[i,3]-1.96*external.f[i,4]^0.5), col = "Red", lty = 2)
    abline(v = (external.f[i,3]+1.96*external.f[i,4]^0.5), col = "Red", lty = 2)
  }
  else {
    hist(Em.f[,i], 
         main = paste0(external.f[i,5], ", ", external.f[i,1], ", ", external.f[i,2]),
         col = "Pink")
    abline(v = external.f[i,3], col = "Red")
    abline(v = (external.f[i,3]-1.96*external.f[i,4]^0.5), col = "Red", lty = 2)
    abline(v = (external.f[i,3]+1.96*external.f[i,4]^0.5), col = "Red", lty = 2)
    
  }
}