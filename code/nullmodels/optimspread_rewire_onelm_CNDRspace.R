# test diffusion model in rewired networks

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/rewire_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/optimfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
log.path <- lapply(Mice,function(X) log(colMeans(X,na.rm = T),base=10))

# load structural networks, actual and rewired

W.orig <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
W.rewired <- readMat(paste(params$opdir,'processed/Wrewire.mat',sep=''))

plt.check.rewire <- function(W.orig,W.rw,ttl){
  p1 <- ggplot() + geom_point(aes(x=rowSums(W.orig),y=rowSums(W.rw)),alpha=0.5,stroke=0) + xlab('Out Degree (Actual)') + ylab('Out Degree (Rewire)') + ggtitle(ttl) + 
    theme_bw()+theme(plot.title = element_text(hjust=0.5),text=element_text(size=6))
  p2 <- ggplot() + geom_point(aes(x=colSums(W.orig),y=colSums(W.rw)),alpha=0.5,stroke=0) + xlab('In Degree (Actual)') + ylab('In Degree (Rewire)') + ggtitle(ttl) + 
    theme_bw()+theme(plot.title = element_text(hjust=0.5),text=element_text(size=6))
  return(plot_grid(plotlist = list(p1,p2),nrow=1))
}
p.in <- plt.check.rewire(W.orig,W.rewired$Wrw.in,'In-Degree Preserved')
p.out <- plt.check.rewire(W.orig,W.rewired$Wrw.out,'Out-Degree Preserved')
p.all <- plot_grid(plotlist= list(p.in,p.out),nrow=2)
ggsave(p.all,filename = paste(savedir,'ConfirmRewiring.pdf',sep=''),
       units = 'cm',height = 9,width = 9,useDingbats=FALSE)

null.list <- list(list(name='InDegreePreserved',W=W.rewired$Wrw.in),list(name='OutDegreePreserved',W=W.rewired$Wrw.out))
for(null.model in null.list){
  null.name <- null.model$name
  W <- null.model$W
  
  L.out.retro <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='retro') # compute out-degreee Laplacian for retrograde connectivity only (not weighted by Snca)
  L.out.antero <- get.Lout(W,rep(1,n.regions.ABA),ant.ret='antero') # compute out-degreee Laplacian for anterograde connectivity only (not weighted by Snca)
  
  #########################################################################################
  ### Use optimization to jointly fit two time scaling parameter on average of all mice ###
  #########################################################################################
  
  Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
  scipy.linalg <- reticulate::import('scipy.linalg') # import scipy matrix exponential function because it's faster
  # c.retro, c.antero: use values from independent fit to initialize
  load(file=paste0(params$opdir,'diffmodel/bidirectional/',paste0(injection.site,collapse='-'),'_independentfit/',grp,'CNDRSpaceIndependentBidirectionalFit_params.RData'))
  params.opt <- c(c.Grp.retro,c.Grp.antero)
  ctrl <- list(fnscale=-1) # maximize the objective function instead of minimizing (default)
  
  params.opt.fit <- optim(params.opt,c.CNDRspace.objective,control = ctrl, method='L-BFGS-B',lower=c(10e-7,10e-7), # optimization. c's must be > 0
                          log.path=log.path,tps=tps,L.out.retro=L.out.retro,L.out.antero=L.out.antero,
                          Xo=Xo,ABA.to.CNDR.key=ABA.to.CNDR.key,fxn =scipy.linalg$expm,one.lm=TRUE) # static inputs
  
  # extract parameters from 
  c.Grp.retro <- params.opt.fit$par[1]
  c.Grp.antero <- params.opt.fit$par[2]
  
  # save data
  Xt.Grp.retro <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.retro,Xo,c.Grp.retro,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
  Xt.Grp.antero <- do.call('cbind',lapply(tps, function(t) log(quiet(map.ABA.to.CNDR(predict.Lout(L.out.antero,Xo,c.Grp.antero,t,fxn=scipy.linalg$expm),path.names,ABA.to.CNDR.key)),base=10))) # predict pathology using connectivity, time constant, and seed
  list[m,e,m.fits,df] <- lm.mask.ant.ret.all(log.path,10^Xt.Grp.retro,10^Xt.Grp.antero) # undo log10 because this function automatically computes it
  save(df,c.Grp.antero,c.Grp.retro,m,m.fits,file = paste(savedir,grp,'CNDRSpace',null.name,'BidirectionalOptim_data.RData',sep=''))
}  
