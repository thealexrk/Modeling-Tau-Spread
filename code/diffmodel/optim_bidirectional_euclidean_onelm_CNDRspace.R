#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_euclidean_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
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

# load distance matrix
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
unit.test(all(region.names == region.names.dmat),'Region names match','ERROR: REGION NAMES DO NOT MATCH')

# exclude missing regions from path and matrix (no data on subiculum)
missing.regions <- unname(which(colSums(is.na(D.mat)) == nrow(D.mat)-1)) # remove regions with all connections missing
region.names.mdl <- region.names.dmat[-missing.regions] # exclude missing regions from model
D.mat <- D.mat[-missing.regions,-missing.regions]
D.mat <- D.mat^-1 # scale so closer regions transmit more path
D.mat[which(diag(nrow(D.mat))==1)] <- 0
L.euclidean <- get.Lout(D.mat,rep(1,nrow(D.mat))) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

L.out.retro <- get.Lout(W[-missing.regions,-missing.regions],rep(1,n.regions.ABA-length(missing.regions)),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W[-missing.regions,-missing.regions],rep(1,n.regions.ABA-length(missing.regions)),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)


#########################################################################################
### Use optimization to jointly fit two time scaling parameter on average of all mice ###
#########################################################################################

Xo <- get.Xo(region.names.mdl,injection.site) # seed pathology in iCPu
injection.site.CNDR <- unname(params$injection.site.CNDR[injection.site]) # convert injection site from ABA to CNDR
scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
#params.opt <- c(0.006070303,0.02223111) # c.retro, c.antero: use values from independent fit to initialize
#params.opt <- c(0.01,0.01) # c.retro, c.antero: use values from independent fit to initialize
load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
params.opt <- c(c.Grp.retro,c.Grp.antero,0.2)
ctrl <- list(fnscale=-1) # minimize objective function

params.opt.fit <- optim(params.opt,c.CNDRspace.euc.objective,control = ctrl, method = 'L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                        log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,L.euclidean=L.euclidean,
                        Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm,one.lm=TRUE,excl.inj.CNDR=injection.site.CNDR) # static inputs

# extract parameters from 
c.Grp.retro <- params.opt.fit$par[1]
c.Grp.antero <- params.opt.fit$par[2]
c.Grp.euclidean <- params.opt.fit$par[3]
# save data
Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
Xt.Grp.euclidean <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.euclidean,Xo,c.Grp.euclidean,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
list[m,e,m.fits,df] <- lm.mask.ant.ret.euc.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero,10^Xt.Grp.euclidean) # undo log10 because this function automatically computes it
save(df,c.Grp.antero,c.Grp.retro,m,m.fits,file = paste(savedir,grp,'CNDRSpaceBidirectionalOptimOneLM_data.RData',sep=''))

# compute vulnerability
df.rnames <- lapply(df,rownames)
vulnerability <- lapply(df.rnames, function(rn) residuals(m)[rn])
pred <- lapply(df.rnames, function(rn) m$fitted.values[rn])
# initialize empty data frame with all CNDR regions and time points
df.vuln <- df.pred <- data.frame(matrix(ncol=length(tps),nrow=n.regions.CNDR, dimnames=list(path.names, paste(tps,'MPI'))),check.names = FALSE)
for(t in 1:length(tps)){
  regions.mpi <- names(vulnerability[[t]]) # get the regions that had non-zero pathology at each tp (so log can be computed)
  regions.mpi.path.names <- substr(regions.mpi,1,nchar(regions.mpi)-2) # get names of regions without "_1" ("_$MPI") appended
  df.vuln[regions.mpi.path.names,paste(tps[t],'MPI')] <- vulnerability[[t]][regions.mpi]
  df.pred[regions.mpi.path.names,paste(tps[t],'MPI')] <- pred[[t]][regions.mpi]
}
tp.excl <- '1 MPI' # exclude 1 MPI - see diffmodel/vuln_hemi_time.R
hemi.average.vuln <- hemi.average(df.vuln[,-which(colnames(df.vuln)==tp.excl)])
save(hemi.average.vuln,df.vuln,df.pred,file = paste0(savedir,grp,'vulnerability_bidirectional_euclidean_hemiaverage_exclude',tp.excl,'.RData'))

