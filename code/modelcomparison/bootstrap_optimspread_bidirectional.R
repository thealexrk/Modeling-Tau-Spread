#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/modelcomparisonfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
# bootstrap log pathology
nboot <- 500
log.path.boot <- bootstrap.path.tps(Mice,nboot)

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

# set up controls for all models
Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
# c.retro, c.antero: use values from bidirectional fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
ctrl.optim <- list(fnscale=-1) # set controls for optim function -- maximize the objective function instead of minimizing (default)
ctrl <- list(Retrograde=L.out.retro,Anterograde=L.out.antero,
             Xo=Xo,fxn=scipy.linalg$expm,ABA.to.CNDR.key=ABA.to.CNDR.key,tps=tps,
             c.r.a.init=c(c.Grp.retro,c.Grp.antero),
             ctrl.optim=ctrl.optim,
             one.lm=FALSE)

# Loop through model types, train on train sets, evaluate on train and test sets
mdl.names <- c('BidirectionalOneLM')
results <- list()
for(mdl.name in mdl.names){
  print(paste('training and testing',mdl.name))
  results[[mdl.name]] <- list(train=NULL)
  for(REP in 1:nboot){
    print(paste('REP',REP))
    X.train <- log.path.boot[[REP]]
    m.out <- spread.model.train(X.train,mdl.name,ctrl)
    results[[mdl.name]]$train[[REP]] <- spread.model.eval(m.out,X.train)
  }
}

# save data
save(results,Mice,file = paste(savedir,grp,'CNDRSpaceModelComparison_Bootstrap.RData',sep=''))
