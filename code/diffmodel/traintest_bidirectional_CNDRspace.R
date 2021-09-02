#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_traintest/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
#params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
#params.opt <- c(0.01,0.01) 
# c.retro, c.antero: use values from bidirectional fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
params.opt <- c(c.Grp.retro,c.Grp.antero)
ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default)

# load path data
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
# mark training samples for each time point
nreps <- 100
tf <- 0.5 # fraction of sample to use for training
train.idx <- lapply(1:nreps, function(R) lapply(Mice, function(X) sample(1:nrow(X),size = tf*nrow(X))))
# split into non-overlapping training and testing sets - get list of lists of mean pathology
log.path.train <- lapply(train.idx, function(t.idx) mapply(function(X,t.idx) {list(log(colMeans(X[t.idx,],na.rm=T),base=10))}, X=Mice,t.idx=t.idx))
log.path.test <- lapply(train.idx, function(t.idx) mapply(function(X,t.idx) {list(log(colMeans(X[-t.idx,],na.rm=T),base=10))}, X=Mice,t.idx=t.idx))

results <- list()
for(REP in 1:nreps){
  X.train <- log.path.train[[REP]] # optimize model on training data
  params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, lower=c(10e-7,10e-7), # optimization. c's must be > 0
                          log.path=X.train,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                          Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm) # static inputs  
  # extract parameters from optimization output
  c.train.retro <- params.opt.fit$par[1]
  c.train.antero <- params.opt.fit$par[2]
  
  # retrieve linear models for antero and retro
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.train.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.train.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  df.train <- lapply(1:length(tps), function(t) data.frame(path = X.train[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
  m.train <- lapply(df.train, function(df.i) lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i))) # fit linear regression at each time point to combine anterograde and retrograde
  train.fits.r <- sapply(m.train, function(m.i) cor(m.i$fitted.values,m.i$model$path)) # extract fits as pearson r
  train.fits.sse <- sapply(m.train, function(m.i) sum(residuals(m.i)^2))
  # apply to test set and measure fit
  X.test <- log.path.test[[REP]]
  df.test <- lapply(1:length(tps), function(t) data.frame(path = X.test[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
  df.test <- lapply(1:length(tps), function(t.) inf.nan.mask(cbind(df.test[[t.]],pred=predict(m.train[[t.]],df.test[[t.]]))) )
  test.fits.r <- sapply(df.test, function(df) cor(df$path,df$pred)) # extract fits as pearson r
  test.fits.sse <- sapply(df.test, function(df) sum((df$path-df$pred)^2)) # extract fits as pearson r
  results[[REP]] <- list(train.idx=train.idx[[REP]], c.train.retro=c.train.retro,c.train.antero=c.train.antero,
                         m.train=m.train,train.fits.r=train.fits.r,train.fits.sse=train.fits.sse,test.fits.r=test.fits.r,test.fits.sse=test.fits.sse)
}

# save data
save(results,Mice,file = paste(savedir,grp,'CNDRSpaceBidirectionalOptim_TrainTest.RData',sep=''))
