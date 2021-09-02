# change directory to example folder
rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
ex.dir <- paste0(basedir,'example/')
setwd(basedir)

# load packages
source('code/misc/packages.R') # see for all of the dependent packages

# load functions
source('code/misc/plottingfxns.R') # I will use some helper functions to provide easy visualization

#######################################################################
# first step: let's develop some intuition of what a linear system is #
#######################################################################

# load the mouse structural connectome. this is an nxn matrix, where 1:(n/2) is 
# one hemisphere, and (n/2):n is the other hemisphere
# This represents the matrix W in the example document

# We are going to implement equation 1: dx/dt = -Lx
# W: connectivity matrix, nxn regions, oriented so that the region in row i sends a projection of strength ij to column j
# D: matrix of zeros whose diagonal contains the row sums (out-degree)
# L: out-degree graph Laplacian of W, which is D - W 
# x: vector specifying how much pathology is in each region, nx1
# t: time. the values of time are arbitrary, because the units of W don't have time in them
# with Allen Brain data, the units of W are some sort of normalized fluorescence intensity, not 1/months

# using L instead of W guarantees the system is stable

load('example/pathdata.RData') # load region names
W <- readMat(paste0(ex.dir,'W.mat'))$W 
n <- nrow(W)
D <- diag(rowSums(W))
L <- D - W

# now let's inject pathology into CA1
x0 <- matrix(0, nrow = n,ncol = 1) # initialize x with all 0's. means no pathology in the brain
injection.site <- c(25) # 27=CA3, 41=DG. can add more injection sites into this vector as in c(25,27,41)
x0[injection.site] <- 1 # now there is 1 unit of pathology in CA1

# now let's simulate this system using equation 7 from the example pdf
# x(t) = e^(-Lt)*x(t=0)

t.n <- 20 # number of time points
t.rng <- seq(0,5,length.out = t.n) # define a range of time in a vector. remember, values of time are arbitrary here
x.t <- matrix(NA, nrow = n, ncol = t.n) # make a matrix to store pathology in whole brain as function of time (x(t))

expm.fxn <- reticulate::import('scipy.linalg')$expm # python's expm function is much faster than R's
# expm.fxn <- expm # if the python function doesn't work, use this line instead, and perhaps delete some nodes in the system

for(t.idx in 1:length(t.rng)){
  print(paste(t.idx,'out of',t.n))
  t. <- t.rng[t.idx]
  x.t[,t.idx] <- expm.fxn(-L*t.)%*%x0 # %*% is R's notation for matrix multiplication. * does elementwise multiplication
}

# now let's see what our system looks like, for the first 10 regions

p <- plot.Xt(x.t,t.rng) + theme(legend.position = 'none') +
  xlab('time (a.u.)') + ylab('Pathology (a.u.)')
p

# you can see the injection site starts at 1 and decays 
# while the other regions peak and decay at different rates

roi.select <- which(rowSums(x.t) > quantile(rowSums(x.t),.95))
p <- plot.Xt(x.t[roi.select,],t.rng) + theme(legend.position = 'none') +
  xlab('time (a.u.)') + ylab('Pathology (a.u.)')
p

# print names of regions that have a relatively strong response to 
# injection of pathology into CA1
print(region.names[roi.select])

#########################################
### Fitting a retrograde spread model ###
#########################################

# But what about this issue of arbitrary time? in a real experiment, we
# measure pathology at specific time points (e.g. 1,3,6 months post inj.)
# There's no reason why the units of W should correspond to our chosen time scale
# To address this potential mismatch, we modify equation 7 by introducing another
# parameter, the time constant c:

# x(t) = e^(-L*c*t)*x(t=0)

# we will use the data to infer what c should be
# our path data is measured in 134 CNDR regions
# the connectivity data exists as 213 ABA regions
# I have written functions that convert to and from these two region definitions
# the fitting has 3 key steps:
# for each time constant value:
# 1. in ABA space, generate predicted pathology values using eq 7 at the relevant time points
# 2. convert that prediction to CNDR space
# 3. compare CNDR space pathology to CNDR prediction
# pick the time constant value that produces the best fit
# the fit is measured as the best spatial correlation between
# predicted and actual, averaged over the measured time points

# This conversion would not be necessary if you measured pathology in ABA space
# or measured connectivity in CNDR space.

# the code below is borrowed from code/diffmodel/analyze_retrogradespread.R

# get mean pathology for each time point
tps <- c(1,3,6,9) 
grp <- 'NTG'
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

# Fit time scaling parameter on average of all mice

# scaling parameter: select range of values. choosing large range will be slow.
# choosing a coarse range will be imprecise. choosing a small range may omit the optimal value

c.rng <- seq(1e-5,0.2,length.out = 50) 
log.path <- lapply(Grp.mean, function(X) log(X,base=10)) # first, we log transfirm the pathology data

# see this function in misc/fitfxns.R. it implements the steps above
# these models don't accurately predict pathology at the injection site
# it forecasts pathology decreases over time when in reality local misfolded leads
# it to increase. So we exclude the injection site from model fitting.
source('code/misc/fitfxns.R')
Xo <- get.Xo(region.names,'iCA1') # seed pathology in CA1
list[c.Grp,Xt.sweep] <- c.CNDRspace.fit(log.path,tps,L,Xo,c.rng,ABA.to.CNDR.key,excl.inj.CNDR = 'iCA1')
print(c.Grp) # this is your optimized time constant value

###############################
### Computing vulenrability ###
###############################

# in our work, we typically define vulnerability as actual - predicted pathology
# positive values indicate a region with more pathology than expected
# negative values indicate less pathology than expected, i.e. resilience
# we will use 2 steps below to get vulnerability
# 1. compute predicted values using optimal time constant
# 2. use univariate linear regression at each time point
# 3. take the residuals of those models to be vulnerability

# predict pathology using connectivity, identified time constant, and seed
Xt.Grp <- do.call('cbind',lapply(tps, function(t) 
  log(quiet(map.ABA.to.CNDR(predict.Lout(L,Xo,c.Grp,t),path.names,ABA.to.CNDR.key)), base = 10))) 
# combine observed pathology and predicted pathology into one data frame at each time point
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred = Xt.Grp[,t,drop=FALSE]))
# fit a linear model at each time point to find deviations between actual and expected
# note that this method gives you a different intercept at each time point, which could be considered a "global" shift in pathology
# specific to each time point
m <- lapply(df, function(df.i) lm(path~pred,data=inf.nan.mask(df.i))) 
# take residuals of those models to get vulnerability at each time point
vulnerability <- lapply(m,residuals)

############################################
### Fitting a bidirectional spread model ###
############################################

# what if pathology spreads in both anterograde and retrograde directions?
# then we can propose a more complex model:
# x(t) = b0 + ba*[e^(-La*ca*t)*x0] + br*[e^(-Lr*cr*t)*x0]

# b0: intercept
# ba: importance of anterograde spread
# br: importance of retrograde spread
# ca: anterograde time constant
# cr: retrograde time constant
# La: anterograde laplacian
# Lr: retrograde laplacian
# we fit this model in a two step process known as alternating optimization

source('code/misc/optimfxns.R')
params.opt <- c(0.006070303,0.02223111)
ctrl <- list(fnscale=-1) # minimize objective function

L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                        log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                        Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =expm.fxn,one.lm=TRUE,excl.inj.CNDR='iCA1') # static inputs

# extract parameters from fit
c.Grp.retro <- params.opt.fit$par[1]
c.Grp.antero <- params.opt.fit$par[2]
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=expm.fxn),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=expm.fxn),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
list[m,e,m.fits,df] <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it

# display summary of model
print(summary(lm.beta(m)))

# vulnerability for all regions at each time point is simply residuals of this model
residuals(m)

