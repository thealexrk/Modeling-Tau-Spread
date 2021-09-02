# This script assesses specificity of model to the actual seed site

rm(list=setdiff(ls(),c('params','grp','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'nullmodels/seedspec_multi/',paste0(injection.site,collapse='-'),'/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/plottingfxns.R')

#################
### Load data ###
#################

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(file = paste(savedir,grp,'AlternateSeedFits_OneLM.RData',sep=''))
load(file = paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',grp,'CNDRSpaceBidirectionalOptimOneLM_data.RData',sep=''))

# compute p-value as probability that other seeds predict observed path better than iCPu seed

null.cors <- seed.fits

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
inj.cor <- m.fits

seed.region.pval <- colMeans(do.call(rbind,lapply(1:nrow(null.cors), function(X) inj.cor)) < null.cors)

#print(paste('Better fit than ',paste0(params$injection.site,collapse=','),':',region.names[which(seed.fits[region.names == 'iCP'] < seed.fits[region.names != 'iCP'])]))

months <- paste(params$tps,'MPI')
months.mat <- do.call(rbind,lapply(1:nrow(null.cors), function(X) months))
p.labs <- paste('p =',signif(seed.region.pval,2))
p.labs[seed.region.pval ==0] <- paste0('p < ',1/nrow(null.cors))
p.null.seeds <- ggplot() + 
  geom_jitter(aes(x=as.vector(months.mat),y = as.vector(null.cors)),color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1, width=0.25) +
  geom_point(aes(x=months,y=as.numeric(inj.cor)),shape = 18,color = 'black',size=2) + 
  geom_text(aes(x=months,y=0.8,label = p.labs),size=2.5) +
  xlab('') + ylab('Fit') + ggtitle('Actual vs. Random Seed') + theme_bw() +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8)) +
  theme(axis.text.x = element_text(size=8)) + theme(axis.text.y = element_text(size=8))
p.null.seeds

ggsave(p.null.seeds,filename=paste(savedir,grp,'SeedSpecificity_OneLM.pdf',sep=''),
	width=7,height=5,units='cm')

# relate alternate fit to connectivity and distance
W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
load(file=paste0(params$opdir,'processed/ABAEuclideanDistanceMatrix.RData'))
dimnames(W) <- list(region.names,region.names)
dimnames(D.mat) <- list(region.names,region.names)

dfun <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
dfun <- function(x1, x2) cor(x1,x2,method='pearson')
inproj.sim <- sapply(alt.seed.sites,function(sites) dfun(rowMeans(W[,sites]),rowMeans(W[,injection.site]))) # correlation with mean connectivity
outproj.sim <- sapply(alt.seed.sites,function(sites) dfun(colMeans(W[sites,]),colMeans(W[injection.site,])))
inproj.sim <- sapply(alt.seed.sites,function(sites) max(dfun(W[,sites],W[,injection.site]))) # mean correlation between connectivity
outproj.sim <- sapply(alt.seed.sites,function(sites) max(dfun(t(W[sites,]),t(W[injection.site,])))) # mean correlation between connectivity
dfun <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
inj.dist <- sapply(alt.seed.sites, function(sites) mean(D.mat[injection.site,sites]))

tps <- params$tps
# fit gams at each time point to compare relationship between null site performance
am <- lapply(1:length(tps), function(t.) gam(null.cors[,t.] ~ s(inproj.sim,k=3) + s(outproj.sim,k=3) + s(inj.dist,k=3),method = 'REML'))
p.list <- lapply(am,function(am.i) p.xy(x=am.i$fitted.values,y=am.i$y,xlab = 'Connectivity Similarity',ylab='Fit to Random Sites',col=wes_palettes$Zissou1[3],alpha=0.5,sm.method = 'gam',formula =y~s(x,k=3)))
p <- plot_grid(plotlist=p.list,align='hv',nrow=1)
ggsave(p,filename=paste(savedir,grp,'RandomSeedFitsVsConnectivity_OneLM.pdf',sep=''),
       height=3.75,width=3.75*length(tps),units='cm')
