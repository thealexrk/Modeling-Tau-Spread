#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
# bootstrap log pathology
log.path.boot <- bootstrap.path.tps(Mice,nboot = 500)

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

# load distance matrix
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
unit.test(all(region.names == region.names.dmat),'Region names match','ERROR: REGION NAMES DO NOT MATCH')

# exclude missing regions from path and matrix (no data on subiculum)
missing.regions <- unname(which(colSums(is.na(D.mat)) == nrow(D.mat)-1)) # remove regions with all connections missing
region.names.mdl <- region.names.dmat[-missing.regions] # exclude missing regions from model
D.mat <- D.mat[-missing.regions,-missing.regions]
D.mat <- D.mat^-1 # scale so closer regions transmit more path
D.mat[which(diag(nrow(D.mat))==1)] <- 0
L.out.D <- get.Lout(D.mat,rep(1,nrow(D.mat))) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

#######################################################################
### Independently fit time constant for antero and retro separately ###
#######################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology in specified injection site
Xo.D <- get.Xo(region.names.mdl,injection.site) # seed pathology in specified injection site using regions with coordinates available
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter

nboot <- length(log.path.boot)
# test these 
mdls.test <- list(Retrograde=list(L=L.out.retro,Xo=Xo),
                  Anterograde=list(L=L.out.antero,Xo=Xo),
                  Euclidean=list(L=L.out.D,Xo=Xo.D))
results <- lapply(mdls.test,function(X) list()) # intialize named results list
df.resid.init <- data.frame(matrix(ncol=length(tps),nrow=n.regions.CNDR, dimnames=list(path.names, paste(tps,'MPI'))),check.names = FALSE) # initialize data frame to hold residuals for every region (even if some are NaN b/c of 0 pathology)
for(mdl.name in names(mdls.test)){
  L.out.mdl <- mdls.test[[mdl.name]]$L # extract connectivity matrix for specified model
  Xo.mdl <- mdls.test[[mdl.name]]$Xo # extract Xo seed site vector for specified model (only differs with euclidean matrix)
  results[[mdl.name]] <- lapply(1:nboot,function(X) list()) # initialize results list for each bootstrap

  for(REP in 1:nboot){
    print(paste('REP',REP))

    X.train <- log.path.boot[[REP]] # optimize time constant for specified model on training data
    list[c.train,Xt.sweep] <- c.CNDRspace.fit(X.train,tps,L.out.mdl,Xo.mdl,c.rng,ABA.to.CNDR.key)  

    # get predicted values
    Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.mdl,Xo.mdl,c.train,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed    
    # fit linear models
    # concatenate pathology with predicted values
    df.train <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred= Xt.Grp[,t,drop=FALSE]))
    m.train <- lapply(df.train, function(df.i) lm(path~pred,data=inf.nan.mask(df.i))) # fit linear regression at each time point based on specified model prediction
    train.fits.mse <- sapply(m.train, function(m.i) mean(residuals(m.i)^2)) # calculated squared deviation from linear model
    train.fits.r <- Xt.sweep[which(c.rng==c.train),] # retrieve correlation from model fit output
    resid <- sapply(m.train, function(m.i) residuals(m.i)) # get residuals of each model
    df.resid <- df.resid.init
    for(tp in 1:length(tps)){df.resid[names(resid[[tp]]),paste(tps[tp],'MPI')] <- resid[[tp]]}
    results[[mdl.name]][[REP]] <- list(c.train=c.train,Xt.Grp=Xt.Grp,resid=df.resid,
                           m.train=m.train,train.fits.r=train.fits.r,train.fits.mse=train.fits.mse)
  }
}

# now construct a new model that is just a linear combo of the independently fit retrograde and anterograde values
mdl.name <- 'BidirectionalIndependent'
results[[mdl.name]] <- list()
for(REP in 1:nboot){
  print(paste('REP',REP))

  X.train <- log.path.boot[[REP]] # load training data for ith rep

  # retrieve predicted values for antero and retro
  Xt.Grp.retro <- results[['Retrograde']][[REP]]$Xt.Grp
  Xt.Grp.antero <- results[['Anterograde']][[REP]]$Xt.Grp
  # concatenate pathology with predicted values
  df.train <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
  m.train <- lapply(df.train, function(df.i) lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i))) # fit linear regression at each time point to combine anterograde and retrograde
  train.fits.mse <- sapply(m.train, function(m.i) mean(residuals(m.i)^2))
  train.fits.r <- sapply(m.train, function(m.i) cor(m.i$fitted.values,m.i$model$path))
  resid <- sapply(m.train, function(m.i) residuals(m.i)) # get residuals of each model
  df.resid <- df.resid.init
  for(tp in 1:length(tps)){df.resid[names(resid[[tp]]),paste(tps[tp],'MPI')] <- resid[[tp]]}

  results[[mdl.name]][[REP]] <- list(c.train=c.train,resid=df.resid,
                         m.train=m.train,train.fits.r=train.fits.r,train.fits.mse=train.fits.mse)
}
mdl.names <- names(results)
save(results,mdl.names,file = paste(savedir,grp,'CNDRSpaceBidirectionalIndependent_Bootstrap.RData',sep=''))
