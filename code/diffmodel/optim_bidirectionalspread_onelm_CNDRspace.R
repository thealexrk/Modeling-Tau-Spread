#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))
log.path <- lapply(Grp.mean, function(X) log(X,base=10))

# load structural networks
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
injection.site.CNDR <- unname(params$injection.site.CNDR[injection.site]) # convert injection site from ABA to CNDR
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
#params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
#params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
params.opt <- c(c.Grp.retro,c.Grp.antero)
ctrl <- list(fnscale=-1) # minimize objective function

params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                        log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                        Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm,one.lm=TRUE,excl.inj.CNDR=NULL) # static inputs

# extract parameters from 
c.Grp.retro <- params.opt.fit$par[1]
c.Grp.antero <- params.opt.fit$par[2]

# save data
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
list[m,e,m.fits,df] <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
save(df,c.Grp.antero,c.Grp.retro,m,m.fits,file = paste(savedir,grp,'CNDRSpaceBidirectionalOptimOneLM_data.RData',sep=''))
