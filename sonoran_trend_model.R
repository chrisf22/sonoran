### Author: Christopher Field
### Date: 5/2022
### Packages required: 'R2jags', 'R2WinBUGS'

### This script specifies a trend model that describes patterns in fire frequency and size in the
### Sonoran Desert that ultimatley will be used to quantify tipping points.

library(R2jags)
library(R2WinBUGS)

# load dataset
fire_dataset <- read.csv("data/fires_sonoran_meters.csv", header=TRUE, stringsAsFactors = FALSE)

# fire frequency
by_year_fun <- function(x) {length(which(fire_dataset$FIRE_YEAR==x))}
fires_per_year <- unlist(lapply(1920:2020, by_year_fun))

plot(fires_per_year, ylim=c(0, 50), pch=16, col=rgb(0, 0, 0, 0.5), bty="n", xlab=c("Year"), ylab=c("Number of fires"), xaxt="n")
axis(side=1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), lab=c(1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020))

# max fire size per year
max_by_year_fun <- function(x) {
  # this approach throws non-missing argument errors but outputs the correct value
  max(max(fire_dataset[fire_dataset$FIRE_YEAR==x, 'SHAPE_Area']), 0)
}
fires_area_year <- unlist(lapply(1920:2020, max_by_year_fun))

# total burned area per year
max_by_year_fun_tot <- function(x) {
  # this approach throws non-missing argument errors but outputs the correct value
  max(sum(fire_dataset[fire_dataset$FIRE_YEAR==x, 'SHAPE_Area']), 0)
}
fires_area_year_total <- unlist(lapply(1920:2020, max_by_year_fun_tot))

plot(fires_area_year, pch=16, col=rgb(0, 0, 0, 0.5), bty="n", xlab=c("Year"), ylab=c("Fire size (units???)"), xaxt="n")
axis(side=1, at=c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), lab=c(1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020))

# get number of years in study period
TT <- length(fires_per_year)

# initial JAGS model to establish a baseline model structure for the number of fires per year
FIRE <- function(){
  int ~ dnorm(0, 0.01)
  # northeast trend
  B ~ dnorm(0, 0.01)
  
  #sd_freq ~ dunif(0, 100)
  #tau_freq <- 1/(sd_freq*sd_freq)
  
  for(t in 1:TT){
    log(lambda[t]) <- int + B*t
    fires_per_year[t] ~ dpois(lambda[t])
  }
}

if (is.R()){
  filename <- file.path(tempdir(), "FIRE.bug")}
write.model(FIRE, filename)
inits <- list(list(int=1, B=0))
data <- list("fires_per_year", "TT")
parameters <- c("int", "B", "lambda")
FIRE <- jags(data=data, inits=inits, parameters.to.save=parameters, filename,
             n.chains=1, n.burnin=50000, n.iter=100000, n.thin=1, DIC=TRUE)

hist(FIRE$BUGSoutput$sims.array[, , "int"])
hist(FIRE$BUGSoutput$sims.array[, , "B"])

# get posterior predictions and caluclate Pearson's and deviance residuals
lambda <- mat.or.vec(TT, 1)
resids_lambda <- mat.or.vec(TT, 1)
dev_resids <- mat.or.vec(TT, 1)
for(t in 1:TT){
  lambda[t] <- mean(FIRE$BUGSoutput$sims.array[, , paste("lambda[", t, "]", sep="")])
  resids_lambda[t] <- (lambda[t] - fires_per_year[t])/sqrt(lambda[t])
  dev_resids[t] <- sign(fires_per_year[t] - lambda[t])*sqrt(2*fires_per_year[t]*max(log(fires_per_year[t]/lambda[t]), 0.001) - (fires_per_year[t] - lambda[t]))
}

plot(dev_resids)
plot(resids_lambda)
