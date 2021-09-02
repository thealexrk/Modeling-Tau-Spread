#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/retrograde/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degreee Laplacian for connectivity only

# Fit time scaling parameter on average of all mice
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))
Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
list[c.Grp,Xt.sweep] <- c.CNDRspace.fit(log.path,tps,L.out,Xo,c.rng,ABA.to.CNDR.key,excl.inj.CNDR = params$excl.inj.CNDR)

Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out,Xo,c.Grp,t),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred = Xt.Grp[,t,drop=FALSE]))
save(df,c.Grp,Xt.sweep,file = paste(savedir,grp,'CNDRSpaceFit_data.RData',sep=''))
