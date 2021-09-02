#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_traintest/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# load results from many split samples
load(file = paste(savedir,grp,'CNDRSpaceBidirectionalOptim_TrainTest.RData',sep=''))

# plot distribution of out-of-sample fits
oos.r <- as.data.frame(do.call(rbind,lapply(results, function(R) R$test.fits.r)))
is.r <- as.data.frame(do.call(rbind,lapply(results, function(R) R$train.fits.r)))
df.plt <- rbind(cbind(oos.r,OI='Test',stringsAsFactors=F),cbind(is.r,OI='Train',stringsAsFactors=F))
colnames(df.plt) <- c(paste(params$tps,'MPI'),'OI')

df.plt <- collapse.columns(df.plt,cnames = paste(params$tps,'MPI'),groupby = 'OI')

p <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25) + theme_classic() +
  ylab('Pearson r') + xlab('') + scale_fill_brewer(palette='Set3',name='') + scale_y_continuous(limits=c(0,NA))+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(p,filename = paste(savedir,grp,'TrainVsTestFits_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 5,useDingbats=FALSE)

# plot distribution of time constants

c.vals <- t(sapply(results,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro)))
df.plt <- collapse.columns(as.data.frame(c.vals))
subset.list <- lapply(c(Anterograde='Anterograde',Retrograde='Retrograde'), function(n) df.plt$values[df.plt$names == n])
subset.stats <- lapply(subset.list, function(X) c(mean=mean(X),sd=sd(X),median=median(X),cv=sd(X)/mean(X)))
sink(paste0(savedir,grp,"TimeConstantStats.txt"))
print(subset.stats)
sink()

p <- ggplot(df.plt) + geom_jitter(aes(x=names,y=values,color=names),size=0.25) + theme_classic() +
  ylab('Diffusivity Constant') + xlab('') + scale_y_continuous(limits=c(0,NA))+scale_color_brewer(palette='Set2',name='')+
  theme(text=element_text(size=8),legend.position='none',axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(p,filename = paste(savedir,grp,'TimeConstants_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 3,useDingbats=FALSE)

# plot distribution of anterograde and retrograde betas

MPI.names <- paste(params$tps,'MPI')
ant.coefs <-as.data.frame(do.call(rbind,lapply(results, function(R) setNames(sapply(R$m.train, function(m.i) coefficients(lm.beta(m.i))['pred.antero']),MPI.names))))
ret.coefs <-as.data.frame(do.call(rbind,lapply(results, function(R) setNames(sapply(R$m.train, function(m.i) coefficients(lm.beta(m.i))['pred.retro']),MPI.names))))
df.plt <- rbind(cbind(ant.coefs,OI='Anterograde',stringsAsFactors=F),cbind(ret.coefs,OI='Retrograde',stringsAsFactors=F))
colnames(df.plt) <- c(MPI.names,'OI')
df.plt <- collapse.columns(df.plt,cnames = MPI.names,groupby = 'OI')
p.lab <- sapply(MPI.names, function(M) pval.2tail.np(0,df.plt$values[df.plt$names==M & df.plt$group=='Anterograde'] - df.plt$values[df.plt$names==M & df.plt$group=='Retrograde']))
p.lab <- sapply(p.lab, function(P) paste0('p = ',signif(P,2)))

p <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size=0.25) + theme_classic() +
  annotate(geom='text',x=MPI.names,y=max(df.plt$values),label=p.lab,size=1)+
  ylab('Standardized Beta') + xlab('') + scale_y_continuous(limits=c(0,NA))+scale_fill_brewer(palette='Set2',name='')+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(p,filename = paste(savedir,grp,'AnterogradeRetrogradeTrainBetas_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 6,useDingbats=FALSE)

