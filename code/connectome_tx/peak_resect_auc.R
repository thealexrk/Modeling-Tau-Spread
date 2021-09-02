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
Xo <- get.Xo(region.names,params$injection.site) # seed pathology in iCPu

L.out <- get.Lout(W,rep(1,n.regions.ABA)) # compute out-degreee Laplacian for connectivity only (not weighted by Snca)
list[c.Grp,Xt.sweep] <- c.CNDRspace.fit(log.path,tps,L.out,Xo,c.rng,ABA.to.CNDR.key)
r.cgrp <- Xt.sweep[which(c.rng==c.Grp),] # overall fit to whole network

## fit model again while deleting 1 connection to injection site at a time

ant.ret <- 'retro'
n.t <- 100 # number of time points to simulate
T. <- 15 # T, end of model time
t.rng <- seq(0.01,T.,length.out = n.t)
c.Grp <- 1
Xt.c <- predict.Lout(L.out,Xo,c.Grp,t=t.rng) # for full model get mode state x(t)
p <- plot.Xt(Xt.conn[[region]],t.rng) + theme(legend.position = 'none',text=element_text(size=8))
ggsave(p,filename = paste(savedir,grp,'Xt_plot.pdf',sep=''),
       units = 'cm',height = 5,width =5 )

Xt.conn <- list() # store model state x(t) for each model with a connection deleted
deleted.conn.names <- setdiff(region.names,params$injection.site)
for(region in deleted.conn.names){
  print(paste0('Deleting ',region,': ',which(region==region.names),'/',n.regions.ABA))
  # for retrograde, delete connection starting at "knockout" region and terminating in injection site
  if(ant.ret == 'retro'){W.del <- W; W.del[which(region.names==region),which(region.names == params$injection.site)] <- 0}
  # for anterograde, delete connection starting at injection site and terminating in "knockout" region
  if(ant.ret == 'antero'){W.del <- W; W.del[which(region.names == params$injection.site),which(region.names==region)] <- 0}
  # compute graph laplacian
  L.out <- get.Lout(W.del,rep(1,n.regions.ABA)) # compute out-degree Laplacian for connectivity only (not weighted by Snca)
  # simulate model
  Xt.conn[[region]] <- predict.Lout(L.out,Xo,c.Grp,t=t.rng)
}

# compute effect on peak time
peak.times.full <- row.Which.Max(Xt.c)
peak.times.del <- lapply(Xt.conn,row.Which.Max)
peak.times.diff <- lapply(peak.times.del,function(X) name(t.rng[peak.times.full] - t.rng[X],region.names)) # compute change in peak time
peak.diff.matrix <- do.call('rbind',peak.times.diff)

peak.times.del <- lapply(peak.times.del,function(X) name(t.rng[X],region.names)) # convert from index to time
peak.times.full <- name(t.rng[peak.times.full],region.names) # convert from index to time


# compute effect on total pathology
dt <- mean(diff(t.rng)) # width of time step
aupc.full <- sum(Xt.c)*dt
aupc.del <- sapply(Xt.conn,sum)*dt
in.deg <- name(colSums(W),region.names)[deleted.conn.names] # get network weighted in degree for non-injection site regions
out.deg <- name(rowSums(W),region.names)[deleted.conn.names] # get network weighted out degree for non-injection site regions

# for retrograde, delete connection starting at "knockout" region and terminating in injection site
if(ant.ret == 'retro'){conn.to.inj <- W[-which(region.names == params$injection.site),which(region.names == params$injection.site)]}
# for anterograde, delete connection starting at injection site and terminating in "knockout" region
if(ant.ret == 'antero'){conn.to.inj <- W[which(region.names == params$injection.site),-which(region.names == params$injection.site)]}

save(peak.times.full,peak.times.del,peak.diff.matrix,c.Grp,aupc.full,aupc.del,out.deg,in.deg,conn.to.inj,file = paste(savedir,grp,'FiberResectionPeakAUC',ant.ret,'.RData',sep=''))

