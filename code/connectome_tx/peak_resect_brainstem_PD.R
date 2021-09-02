# how does deleting connections affect overall pathology

#################
### Load data ###
#################

grp <- 'NTG'

rm(list=setdiff(ls(),c('params','grp')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'connectome_tx/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

## get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

## load connectome
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W

################
### Analysis ###
################

## fit full model

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter -- shortening to increase speed
log.path <- lapply(Grp.mean, function(X) log(X,base=10))

# load time constant
load(file=Sys.glob(paste0(params$opdir,'diffmodel/retrograde/*/',grp,'CNDRSpaceFit_data.RData')))
L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)

## fit model again while deleting 1 connection to injection site at a time

ant.ret <- 'retro'
n.t <- 100 # number of time points to simulate
T. <- 15 # T, end of model time
t.rng <- seq(0.01,T.,length.out = n.t)
Xo <- get.Xo(region.names,params$injection.site)
c.Grp <- 1
Xt.c <- predict.Lout(L.out,Xo,c.Grp,t=t.rng) # for full model get mode state x(t)
p <- plot.Xt(Xt.c,t.rng) + theme(legend.position = 'none',text=element_text(size=8))
ggsave(p,filename = paste(savedir,grp,'Xt_plot.pdf',sep=''),
       units = 'cm',height = 5,width =5 )

deleted.sumstats <- peak.times.full <- Xt.c <- list() # store summary statistics for deleting each connection one at a time
brainstem.PD.sites <- list(LC='iPCG',CP='iCP',Amygdala=c('iBLA','iBMA','iLA','iCEA','iPA','iMEA','iAAA','iA'),
                      SN=c('iSNc','iSNr'),Insula=c('iAIp','iAIv','iAId'),NTS='iNTS',XII='iXII')
for(injection.site.name in names(brainstem.PD.sites)){
  injection.site <- brainstem.PD.sites[[injection.site.name]]
  Xt.conn <- list() # store model state x(t) for each model with a connection deleted
  deleted.conn.names <- setdiff(region.names,injection.site)
  Xo <- get.Xo(region.names,injection.site) # seed pathology in iCPu
  Xt.c[[injection.site.name]] <- predict.Lout(L.out,Xo,c.Grp,t=t.rng) # for full model get mode state x(t)
  peak.times.full[[injection.site.name]] <- row.Which.Max(Xt.c[[injection.site.name]]) # find when each region peaks
  for(region in deleted.conn.names){
    print(paste0('Deleting ',region,': ',which(region==region.names),'/',n.regions.ABA))
    # for retrograde, delete connection starting at "knockout" region and terminating in injection site
    if(ant.ret == 'retro'){W.del <- W; W.del[which(region.names==region),which(region.names %in% injection.site)] <- 0}
    # for anterograde, delete connection starting at injection site and terminating in "knockout" region
    if(ant.ret == 'antero'){W.del <- W; W.del[which(region.names %in% injection.site),which(region.names==region)] <- 0}
    # compute graph laplacian
    L.out <- get.Lout(W.del,rep(1,n.regions.ABA)) # compute out-degree Laplacian for connectivity only (not weighted by Snca)
    # simulate model
    Xt.conn[[injection.site.name]][[region]] <- predict.Lout(L.out,Xo,c.Grp,t=t.rng)
  }

  # compute effect on peak time
  peak.times.del <- lapply(Xt.conn[[injection.site.name]],row.Which.Max)
  peak.times.diff <- lapply(peak.times.del,function(X) name(t.rng[peak.times.full[[injection.site.name]]] - t.rng[X],region.names)) # compute change in peak time
  peak.diff.matrix <- do.call('rbind',peak.times.diff)
  
  peak.times.del <- lapply(peak.times.del,function(X) name(t.rng[X],region.names)) # convert from index to time
  
  # compute effect on total pathology
  dt <- mean(diff(t.rng)) # width of time step
  aupc.del <- sapply(Xt.conn[[injection.site.name]],sum)*dt
  

  # store these summary stats for each injection site
  deleted.sumstats[[injection.site.name]]$peak.diff.matrix <- peak.diff.matrix
  deleted.sumstats[[injection.site.name]]$peak.times.del <- peak.times.del
  deleted.sumstats[[injection.site.name]]$aupc.del <- aupc.del
  # for retrograde, delete connection starting at "knockout" region and terminating in injection site
  if(ant.ret == 'retro'){deleted.sumstats[[injection.site.name]]$conn.to.inj <- W[-which(region.names %in% injection.site),which(region.names %in% injection.site)]}
  # for anterograde, delete connection starting at injection site and terminating in "knockout" region
  if(ant.ret == 'antero'){deleted.sumstats[[injection.site.name]]$conn.to.inj <- W[which(region.names %in% injection.site),-which(region.names %in% injection.site)]}
  names(deleted.sumstats[[injection.site]]$conn.to.inj) <- region.names[which(!region.names %in% injection.site)]

}
in.deg <- name(colSums(W),region.names)[deleted.conn.names] # get network weighted in degree for non-injection site regions
out.deg <- name(rowSums(W),region.names)[deleted.conn.names] # get network weighted out degree for non-injection site regions


n.t <- 100 # number of time points to simulate
T. <- 15 # T, end of model time
t.rng <- seq(0.01,T.,length.out = n.t)
dt <- mean(diff(t.rng))
c.Grp <- 1
L.out <- get.Lout(W,rep(1,n.regions.ABA))
aupc.baseline <- list()
for(injection.site.name in names(brainstem.PD.sites)){
  print(injection.site.name)
  Xo <- get.Xo(region.names,brainstem.PD.sites[[injection.site.name]])
  aupc.baseline[[injection.site.name]] <- sum(predict.Lout(L.out,Xo,c.Grp,t=t.rng))*dt # for full model get mode state x(t)
}

save(deleted.sumstats,out.deg,in.deg,brainstem.PD.sites,aupc.baseline,file = paste(savedir,grp,'FiberResectionPeakAUC',ant.ret,'.RData',sep=''))

