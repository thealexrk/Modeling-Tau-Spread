# This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/seedspec/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')

#################
### Load data ###
#################

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for only GBA time point, 1 month
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
log.path <- lapply(Mice,function(X) log(colMeans(X,na.rm = T),base=10))

# Fit time scaling parameter on average of all mice
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degree Laplacian for connectivity only (not weighted by Snca)

c.rng <- seq(params$c.min,params$c.max*10,length.out = params$c.n*2) # scaling parameter
seed.fits <- matrix(NA,nrow=n.regions.ABA,ncol=length(tps),dimnames=list(region.names,paste(tps,'MPI'))) # preallocate vector to store fits from each seed region
c.altseed <- matrix(NA,nrow=n.regions.ABA,dimnames=list(region.names,NULL))

for(S in 1:n.regions.ABA){
  seed.name <- region.names[S]
  print(paste0('Seed: ',seed.name,', ',S,' out of ',n.regions.ABA))
  Xo <- get.Xo(region.names,seed.name) # seed pathology in each region one at a time
  list[c.Grp,Xt.sweep] <- c.CNDRspace.fit(log.path,tps,L.out,Xo,c.rng,ABA.to.CNDR.key) # fit time constant
  seed.fits[S,] <- Xt.sweep[which(c.rng == c.Grp),]
  c.altseed[S] <- c.Grp
}

save(seed.fits,c.altseed, file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))
