#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/traintest/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/modelcomparisonfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

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

# load path data and get training/testing sets
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
nreps <- 500
tf <- 0.5
list[log.path.train,log.path.test,train.idx] <- train.test.logmeanpath(Mice,nreps = nreps,tf = 0.5)

# set up controls for all models
Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
injection.site.CNDR <- unname(params$injection.site.CNDR[injection.site]) # convert injection site from ABA to CNDR
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
  # c.retro, c.antero: use values from bidirectional fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
ctrl.optim <- list(fnscale=-1) # set controls for optim function -- maximize the objective function instead of minimizing (default)
ctrl.base <- list(Retrograde=L.out.retro,Anterograde=L.out.antero,Euclidean=L.out.D,
             Xo=Xo,fxn=scipy.linalg$expm,ABA.to.CNDR.key=ABA.to.CNDR.key,tps=tps,
             c.r.a.init=c(c.Grp.retro,c.Grp.antero),
             ctrl.optim=ctrl.optim,
             c.rng=seq(params$c.min,params$c.max,length.out = params$c.n),
             one.lm=FALSE)

# Loop through model types, train on train sets, evaluate on train and test sets
mdl.names <- c('Euclidean','Retrograde','Anterograde','Bidirectional','BidirectionalOneLM','BidirectionalOneLMEuclidean')
results <- list()

for(mdl.name in mdl.names){
  ctrl <- ctrl.base
  print(paste('training and testing',mdl.name))
  
  if(grepl('Euclidean',mdl.name)){
    ctrl$Xo <- get.Xo(region.names.mdl,injection.site) # exclude regions with missing coordinate data
    ctrl$excl.inj <- injection.site.CNDR # exclude injection sites which end up being major outliers in this model
    ctrl$c.rng <- ctrl$c.rng + 0.1 # shift c.rng to capture different range of optimal time constants in Euclidean distance model
    if(grepl('Bidirectional',mdl.name)){ # bidirectional euclidean distance model; delete regions with missing coordinates
      ctrl$Retrograde <- get.Lout(W[-missing.regions,-missing.regions],rep(1,length(region.names.mdl)),ant.ret = 'retro')
      ctrl$Anterograde <- get.Lout(W[-missing.regions,-missing.regions],rep(1,length(region.names.mdl)),ant.ret = 'antero')
      ctrl$c.r.a.init <- c(ctrl$c.r.a.init,0.2) # initialize euclidean distance time constant at 0.2
    }
  }
  
  results[[mdl.name]] <- list(train=NULL,test=NULL)
  for(REP in 1:nreps){ 
    print(paste('REP',REP))
    X.train <- log.path.train[[REP]]
    X.test <- log.path.test[[REP]]
    m.out <- spread.model.train(X.train,mdl.name,ctrl)
    results[[mdl.name]]$train[[REP]] <- spread.model.eval(m.out,X.train)
    results[[mdl.name]]$test[[REP]] <- spread.model.eval(m.out,X.test)
  }
}

# save data
save(results,Mice,file = paste(savedir,grp,'CNDRSpaceModelComparison_TrainTest.RData',sep=''))
