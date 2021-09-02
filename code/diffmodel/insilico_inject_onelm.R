#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'insilico_injections_onelm/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# load time constant and model objects for bidirectional, by default using all specified injection sites
datadir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(params$injection.site,collapse='-'),'/',sep='')
load(file=paste0(datadir,grp,'CNDRSpaceBidirectionalOptimOneLM_data.RData'))
tps <- params$tps
# load connectivity matrix, compute antero and retro laplacians
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)

Xo <- get.Xo(region.names,injection.site) # seed pathology in specified in silico site

# get predicted pathology in CNDR space
Xt.Grp.antero <- log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,params$tps),path.names,ABA.to.CNDR.key)), base = 10) # predict pathology using connectivity, time constant, and seed
Xt.Grp.retro <- log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,params$tps),path.names,ABA.to.CNDR.key)), base = 10) # predict pathology using connectivity, time constant, and seed
df <- lapply(1:length(params$tps), function(t) data.frame(x1 = Xt.Grp.retro[,t,drop=FALSE], x2 = Xt.Grp.antero[,t,drop=FALSE]))

insilico.pred <- matrix(NA,nrow=n.regions.CNDR,ncol=length(params$tps),dimnames = list(path.names,paste(tps,'MPI')))

for(t. in 1:length(tps)){df[[t.]][,'pred.retro'] <- 0}
m.pred <- lapply(1:length(tps), function(t.) predict(m,df[[t.]]))
for(t. in 1:length(tps)){ insilico.pred[names(m.pred[[t.]]),t.] <- m.pred[[t.]]}
write.csv(x=insilico.pred,file = paste0(savedir,grp,'InsilicoInjection_',injection.site,'.csv'))

source('code/misc/plottingfxns.R')
#Xt <- do.call(cbind,lapply(1:length(tps), function(t.) scale(insilico.pred[,t.,drop=FALSE],center=T)))

plot.Xt(insilico.pred,tps) +theme(legend.position = 'none')

