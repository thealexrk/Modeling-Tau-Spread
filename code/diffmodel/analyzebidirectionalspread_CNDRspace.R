#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))
Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
list[c.Grp.retro,Xt.sweep.retro] <- c.CNDRspace.fit(log.path,tps,L.out.retro,Xo,c.rng,ABA.to.CNDR.key)
list[c.Grp.antero,Xt.sweep.antero] <- c.CNDRspace.fit(log.path,tps,L.out.antero,Xo,c.rng,ABA.to.CNDR.key)

scipy.linalg <- reticulate::import('scipy.linalg')
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
save(df,c.Grp.antero,c.Grp.retro,Xt.sweep.retro,Xt.sweep.antero,file = paste(savedir,grp,'CNDRSpaceBidirectionalFit_data.RData',sep=''))
save(c.Grp.antero,c.Grp.retro,file = paste(savedir,grp,'CNDRSpaceIndependentBidirectionalFit_params.RData',sep=''))
