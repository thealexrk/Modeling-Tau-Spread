#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/euclidean/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

# load distance matrix
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
unit.test(all(region.names == region.names.dmat),'Region names match','ERROR: REGION NAMES DO NOT MATCH')

# exclude missing regions from path and matrix (no data on subiculum)
missing.regions <- unname(which(colSums(is.na(D.mat)) == nrow(D.mat)-1)) # remove regions with all connections missing
region.names.mdl <- region.names.dmat[-missing.regions] # exclude missing regions from model
D.mat <- D.mat[-missing.regions,-missing.regions]
D.mat <- D.mat^-1 # scale so closer regions transmit more path
D.mat[which(diag(nrow(D.mat))==1)] <- 0 # remove diagonal
L.out <- get.Lout(D.mat,rep(1,nrow(D.mat))) # compute out-degreee Laplacian for connectivity only (not weighted by expression)

# Fit time scaling parameter on average of all mice
#c.rng <- seq(params$c.min/100,params$c.max/100,length.out = params$c.n) # scaling parameter
# adding 0.1 because optimal c is around 0.26 if you exclude injection site
c.rng <- seq(params$c.min+0.1,params$c.max+0.1,length.out = params$c.n) # scaling parameter
c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
log.path <- lapply(Grp.mean, function(X) log(X,base=10))

Xo <- get.Xo(region.names.mdl,injection.site) # seed pathology in iCPu

# c.Grp <- 0.014
# t.rng <- seq(0,15,length.out = 100)
# Xt.c <- predict.Lout(L.out,Xo,c.Grp,t=t.rng) # for full model get mode state x(t)
# source('code/misc/plottingfxns.R')
# p <- plot.Xt(Xt.c,t.rng) + theme(legend.position = 'none',text=element_text(size=8))
# p

injection.site.CNDR <- unname(params$injection.site.CNDR[injection.site]) # convert injection site from ABA to CNDR
list[c.Grp,Xt.sweep] <- c.CNDRspace.fit(log.path,tps,L.out,Xo,c.rng,ABA.to.CNDR.key,excl.inj = injection.site.CNDR)

Xt.Grp <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out,Xo,c.Grp,t),path.names,ABA.to.CNDR.key)), base = 10))) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(tps), function(t) data.frame(path = log.path[[t]], pred = Xt.Grp[,t,drop=FALSE]))
save(df,c.Grp,Xt.sweep,file = paste(savedir,grp,'CNDRSpaceEuclideanDistanceFit_data.RData',sep=''))
