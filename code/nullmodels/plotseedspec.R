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
load(file = paste(savedir,grp,'AlternateSeedFits.RData',sep=''))
load(file = paste(params$opdir,'diffmodel/',grp,'CNDRSpaceFit_data.RData',sep=''))

# compute p-value as probability that other seeds predict observed path better than iCPu seed

null.cors <- seed.fits[!region.names %in% params$injection.site,]
inj.cor <- seed.fits[region.names %in% params$injection.site,]

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n) # scaling parameter
inj.cor <- name(Xt.sweep[which(c.rng ==c.Grp),],paste(params$tps,'MPI'))

seed.region.pval <- colMeans(do.call(rbind,lapply(1:nrow(null.cors), function(X) inj.cor)) < null.cors)

#print(paste('Better fit than ',paste0(params$injection.site,collapse=','),':',region.names[which(seed.fits[region.names == 'iCP'] < seed.fits[region.names != 'iCP'])]))

months <- paste(params$tps,'MPI')
months.mat <- do.call(rbind,lapply(1:nrow(null.cors), function(X) months))
p.null.seeds <- ggplot() + 
  geom_jitter(aes(x=as.vector(months.mat),y = as.vector(null.cors)),color ='#5F4B8B',alpha = 0.5,stroke = 0,size = 1, width=0.25) +
  geom_point(aes(x=months,y=as.numeric(inj.cor)),shape = 18,color = 'black',size=2) + 
  geom_text(aes(x=months,y=0.8,label = paste('p =',signif(seed.region.pval,2))),size=2.5) +
  xlab('') + ylab('Fit') + ggtitle('Actual vs. Random Seed') + theme_bw() +
  theme(text = element_text(size=8),plot.title = element_text(hjust=0.5,size=8)) +
  theme(axis.text.x = element_text(size=8)) + theme(axis.text.y = element_text(size=8))
p.null.seeds

ggsave(p.null.seeds,filename=paste(savedir,grp,paste0(params$injection.site,collapse=','),'SeedSpecificity.pdf',sep=''),
	width=7,height=5,units='cm')
