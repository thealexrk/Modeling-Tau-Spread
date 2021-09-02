# This script assesses specificity of model to the actual seed site
# accounting for the fact that we are using multiple seed sites
# and there may be a benefit of picking spatially contiguous seed sites

rm(list=setdiff(ls(),c('params','grp','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'nullmodels/seedspec_multi/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')

#################
### Load data ###
#################

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for only GBA time point, 1 month
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
log.path <- lapply(Mice,function(X) log(colMeans(X,na.rm = T),base=10))

# load distance matrix
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
unit.test(all(region.names == region.names.dmat),'Region names match','ERROR: REGION NAMES DO NOT MATCH')

nperms <- 500 # test 500 combinations of different seed sites

# identify alternate combinations of seed sites that are as close together as injection.site sites
# original distance volume has pixel size of 25 microns, but here we are just using array indices which are proportional to that
D.mat <- D.mat*(25/1000) # convert distance into mm
dimnames(D.mat) <- list(region.names.dmat,region.names.dmat)
n.sites <- length(injection.site)
d.obs <- sum(tril(D.mat[injection.site,injection.site]))/(n.sites*(n.sites-1)/2) # observed average distance between injection site is sum of submatrix (*1/2 but not needed b/c it's symmetric)

non.seed.sites <- region.names[which(!region.names %in% injection.site) ] # get names of regions that are not the tested injection sites
D.mat.test <- abs(D.mat - d.obs) # Normalize D.mat to the average distance between injection sites

d.thresh <- d.obs*0.1 # only add sites that will change mean distance by 10% of target
alt.seed.sites <- list()
for(P in 1:ceiling(nperms*1.1)){ # do 10% more in case some iterations fail
  # randomly select a region
  alt.seeds <- sample(non.seed.sites,1)
  total.distance <- 0 # keep running total of inter-region distance
  for(k in 2:n.sites){
    new.region.potential <- sort(abs((total.distance+colSums(D.mat[alt.seeds,,drop=FALSE]))/(k*(k-1)/2) - d.obs)) # sort all other regions by average distance to seed site(s)
    new.region.potential <- names(new.region.potential[new.region.potential < d.thresh ])
    # add regions randomly as long as they don't increase total distance above threshold
    if(length(new.region.potential) >0){alt.seeds <- c(alt.seeds,sample(new.region.potential,1))}
    total.distance <- sum(tril(D.mat[alt.seeds,alt.seeds,drop=FALSE]))
  }
  alt.seed.sites[[P]] <- alt.seeds
}
alt.seed.sites <- alt.seed.sites[sapply(alt.seed.sites,length) ==5] # throw away any combos that were too far
alt.seed.sites <- alt.seed.sites[1:nperms] # truncate to specified number of permutations
alt.seed.distances <- sapply(alt.seed.sites, function(X) sum(tril(D.mat[X,X]))/(n.sites*(n.sites-1)/2) )

p <- ggplot() + geom_histogram(aes(x=alt.seed.distances),fill=wes_palettes$Royal2[3]) + geom_vline(xintercept = c(d.obs-d.thresh,d.obs,d.obs+d.thresh),linetype='dashed')+
  theme_bw() + theme(text=element_text(size=8)) + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(limits=c(0,(d.obs+d.thresh)*1.1))+
  xlab('Mean Inter-Seed\nDistance (mm)')
ggsave(p,filename = paste(savedir,'AlternateSeedSiteDistances.pdf',sep=''),
       units = 'cm',height = 4,width = 4,useDingbats=FALSE)

seed.fits <- matrix(NA,nrow=nperms,ncol=length(tps),dimnames=list(NULL,paste(tps,'MPI'))) # preallocate vector to store fits from each seed region
c.retro.altseed <- matrix(NA,nrow=nperms,dimnames=list(NULL,NULL))
c.antero.altseed <- matrix(NA,nrow=nperms,dimnames=list(NULL,NULL))

# prepare connectivity matrices to fit alternative seed sites with bidirectional model
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default)

for(S in 132:nperms){
  seed.names <- alt.seed.sites[[S]]
  print(paste0('Seed: ',paste0(seed.names,collapse='-'),', ',S,' out of ',nperms))
  Xo <- get.Xo(region.names,seed.names) # seed pathology in each set of seed regions
  params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, lower=c(10e-7,10e-7), # optimization. c's must be > 0
                          log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                          Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm) # static inputs
  c.retro.altseed[S] <- params.opt.fit$par[1]
  c.antero.altseed[S] <- params.opt.fit$par[2]
  
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.retro.altseed[S],t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.antero.altseed[S],t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
  df <- lapply(1:length(params$tps), function(t) data.frame(path = log.path[[t]], pred.retro = Xt.Grp.retro[,t,drop=FALSE], pred.antero = Xt.Grp.antero[,t,drop=FALSE]))
  m <- lapply(df, function(df.i) lm(path~pred.retro+pred.antero,data=inf.nan.mask(df.i)))
  seed.fits[S,] <- sapply(m, function(m.i) cor(m.i$fitted.values,m.i$model$path))
}

save(seed.fits,c.retro.altseed,c.antero.altseed,alt.seed.distances,alt.seed.sites, file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))
