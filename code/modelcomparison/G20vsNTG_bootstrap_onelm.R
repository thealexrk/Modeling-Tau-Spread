#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'modelcomparison/bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

###################################################
### load bootstrap data for bidirectional model ###
###################################################

mdl.name <- 'BidirectionalOneLM'
load(file = paste(savedir,'NTGCNDRSpaceModelComparison_Bootstrap.RData',sep=''))
results.NTG <- results[[mdl.name]]$train

load(file = paste(savedir,'G20CNDRSpaceModelComparison_Bootstrap.RData',sep=''))
results.G20 <- results[[mdl.name]]$train
rm(results)

MPI.names <- paste(params$tps,'MPI')
group.colors <- getGroupColors(params$grps)

# 1. compare model fits between g20 and NTG
cfg <- list(fits.mse='MSE',fits.r='Pearson r') # loop through pearson r and MSE
p.fit <- list()
for(f.met in names(cfg)){ # compare fits for MSE and pearson r
  NTG.r <- as.data.frame(do.call(rbind,lapply(results.NTG, function(R) R[[f.met]])))
  G20.r <- as.data.frame(do.call(rbind,lapply(results.G20, function(R) R[[f.met]])))
  df.all <- rbind(cbind(NTG.r,Group='NTG',stringsAsFactors=F),cbind(G20.r,Group='G20',stringsAsFactors=F))
  colnames(df.all) <- c(MPI.names,'Group')
  df.plt <- collapse.columns(df.all,cnames = MPI.names,groupby = 'Group')
  df.plt$group <- factor(df.plt$group,levels = c(params$grps),ordered=T)
  p.lab <- sapply(1:ncol(NTG.r), function(tp) pval.2tail.np(0,NTG.r[,tp]-G20.r[,tp]))
  p.fit[[cfg[[f.met]]]] <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size=0.25) + theme_classic() +
    annotate(geom='text',x=MPI.names,y=Inf,label=ifelse(p.lab<0.05,yes='*',no=''),vjust=1,size=2)+
    ylab(cfg[[f.met]]) + xlab('') + scale_fill_manual(limits=params$grps,values =group.colors,name='') + scale_y_continuous(limits=c(0,NA))+
    theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
  ggsave(p.fit[[cfg[[f.met]]]],filename = paste(savedir,'NTGvsG20Bootstrap',cfg[[f.met]],'_CNDRSpace.pdf',sep=''),
         units = 'cm',height = 4,width = 5,useDingbats=FALSE)
}

# 3. compare time constants between G20 and NTG

NTG.c <- as.data.frame(t(sapply(results.NTG,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
G20.c <- as.data.frame(t(sapply(results.G20,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
df.all <- rbind(cbind(NTG.c,Group='NTG',stringsAsFactors=F),cbind(G20.c,Group='G20',stringsAsFactors=F))
df.plt <- collapse.columns(df.all,cnames=c('Anterograde','Retrograde'),groupby = 'Group')
for(grp in params$grps){
  subset.list <- lapply(c(Anterograde='Anterograde',Retrograde='Retrograde'), function(n) df.plt$values[df.plt$names == n & df.plt$group==grp])
  subset.stats <- lapply(subset.list, function(X) c(mean=mean(X),quantile(X,c(0.025,.975)),median=median(X),cv=sd(X)/mean(X)))
  subset.stats$AnterogradeVsRetrogradeTimeConstant <- pval.np.pub(pval.2tail.np(0,subset.list$Retrograde-subset.list$Anterograde),length(subset.list$Retrograde))
  sink(paste0(savedir,grp,"TimeConstantStats.txt"))
  print(subset.stats)
  sink()
}
p.lab <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0,df.plt$values[df.plt$names==d & df.plt$group == 'NTG']-df.plt$values[df.plt$names==d & df.plt$group == 'G20']))
p.lab <- paste0('p = ',signif(p.lab,2))
df.plt$group <- factor(df.plt$group,levels = c(params$grps),ordered=T)
p.c <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size = 0.25) + theme_classic() +
  annotate(geom='text',x=c('Anterograde','Retrograde'),y=Inf,label=p.lab,vjust=1,size=2)+
  ylab('Diffusion Rate Constant') + xlab('') + scale_y_continuous(limits=c(0,NA))+scale_fill_manual(limits=params$grps,values =group.colors,name='')+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
p.c
ggsave(p.c,filename = paste(savedir,'NTGvsG20TimeConstants_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 5,useDingbats=FALSE)

# 3. compare anterograde and retrograde betas between groups
cfg <- list(NTG=results.NTG,G20=results.G20) # loop through group results
for(grp in params$grps){
  # see code/misc/optimfxns.R, function lm.mask.ant.ret.all: x1 = retrograde, x2 = anterograde
  # ^ -- probably should change it so x1 is just called retro to begin with
  # extract coefficients weighting the importance of anterograde and retrograde betas
  ant.coefs <- as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) R$m.train.coefs['x2','Estimate'])))
  ret.coefs <-as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) R$m.train.coefs['x1','Estimate'])))
  
  df.plt <- rbind(cbind(ant.coefs,OI='Anterograde',stringsAsFactors=F),cbind(ret.coefs,OI='Retrograde',stringsAsFactors=F))
  colnames(df.plt) <- c('Beta','OI')
  df.plt <- collapse.columns(df.plt,cnames = 'Beta',groupby = 'OI')
  df.plt$Group <- grp
  p.avr <- pval.2tail.np(0,df.plt$values[df.plt$group=='Anterograde'] - df.plt$values[df.plt$group=='Retrograde'])
  p.avr <- p.avr
  p.lab.avr <- paste0('p = ',signif(p.avr,2))
  p.lab.avr[p.avr ==0] <- paste0('p < ',signif(1/length(cfg[[grp]]),2)) # don't say p = 0, say < 1/nboots
  cfg[[grp]] <- list(df=df.plt,p.avr=p.lab.avr,grp=grp) # store plot ready df in looped list
}

df.plt <- do.call(rbind,lapply(cfg,function(X) X$df)) # vertically concatenate df across G20 and NTG then plot
# p-value comparing ant vs. ret for each group
p.avr.save <- do.call(rbind,lapply(cfg, function(X) c(X$p.avr,Group=X$grp,Pcrit=0.05/length(X$p.avr))))
write.csv(x=p.avr.save,file=paste0(savedir,'NTGandG20_AnteroVsRetroBetas_Stats.csv'),row.names = F)
# compare the two groups to each other
p.G20.v.NTG <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0, # loop through direction
                  cfg[['NTG']]$df$values[cfg[['NTG']]$df$group==d] - # compare group diff to 0
                    cfg[['G20']]$df$values[cfg[['G20']]$df$group==d]))
p.G20.v.NTG.lab <- data.frame(values=paste('p =',signif(p.G20.v.NTG,2)))
p.G20.v.NTG.lab$group <- c('Anterograde','Retrograde')
df.plt$Group <- factor(df.plt$Group,levels = params$grps,ordered=T)
p <- ggplot() + geom_boxplot(data=df.plt,aes(x=group,y=values,fill=Group),size=0.25,outlier.size=0.25) + 
  theme_classic() +
  geom_text(data=p.G20.v.NTG.lab,aes(x=group,label=values,y=Inf),vjust=1,size=1.8)+
  ylab('Standardized Beta') + xlab('') +scale_fill_manual(limits=params$grps,values =group.colors,name='')+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(p,filename = paste(savedir,'NTGvsG20AnterogradeRetrogradeBetas_Boxplot_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 6,width = 9,useDingbats=FALSE)
write.csv(x=p.G20.v.NTG.lab,file=paste0(savedir,'NTGvsG20_AnteroAndRetroBetas_Stats.csv'),row.names = F)

# save the mean and 95% CI of betas and time constants within each group
NTG.c <- as.data.frame(t(sapply(results.NTG,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
G20.c <- as.data.frame(t(sapply(results.G20,function(R) c(Anterograde=R$c.train.antero,Retrograde=R$c.train.retro))))
df.tc <- rbind(cbind(NTG.c,Group='NTG',stringsAsFactors=F),cbind(G20.c,Group='G20',stringsAsFactors=F))
df.tc <- collapse.columns(df.all,cnames=c('Anterograde','Retrograde'),groupby = 'Group')

beta.tc.stats <- data.frame()

for(a.r in c('Anterograde','Retrograde')){
  for(grp in params$grps){
    mask.beta <- df.plt$group==a.r &df.plt$Group==grp
    beta.tc.stats[paste(grp,a.r,'Beta'),'2.5%'] <- quantile(df.plt$values[mask.beta],0.025)
    beta.tc.stats[paste(grp,a.r,'Beta'),'Median'] <- median(df.plt$values[mask.beta])
    beta.tc.stats[paste(grp,a.r,'Beta'),'97.5%'] <- quantile(df.plt$values[mask.beta],0.975)
  }
}
for(a.r in c('Anterograde','Retrograde')){
  for(grp in params$grps){
    mask.tc <- df.tc$names == a.r & df.tc$group==grp
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'2.5%'] <- quantile(df.tc$values[mask.tc],0.025)
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'Median'] <- median(df.tc$values[mask.tc])
    beta.tc.stats[paste(grp,a.r,'Diffusion Time Constant'),'97.5%'] <- quantile(df.tc$values[mask.tc],0.975)
  }
}
write.csv(x=signif(beta.tc.stats,3),file=paste0(savedir,'NTGvsG20BetaTimeConstantDistributionTable.csv'))

# 4. bootstrap NTG vulnerability values -- not done

NTG.v <- lapply(results.NTG,function(R) R$resid)
# first bootstrap at each time point
NTG.v.by.tp <- lapply(MPI.names, function(tp)
  do.call(cbind,lapply(NTG.v, function(X) X[,tp,drop=F])))
p.v0.by.tp <- lapply(NTG.v.by.tp, function(M.v) # loop through each months vulnerability
  # computed 2-tailed p-value asking if each region's vulnerability differs from 0
  2*apply(cbind(rowMeans(M.v>0,na.rm = T),rowMeans(M.v<0,na.rm=T)),1,min))
p.v0.by.tp <- do.call(cbind,p.v0.by.tp) # concatenate horizontally
# next bootstrap hemisphere average vulnerability
NTG.v.hemi <- do.call(cbind,lapply(results.NTG,function(R) hemi.average(R$resid[,MPI.names[-1]])))
p.v0.hemi <- 2*apply(cbind(rowMeans(NTG.v.hemi>0,na.rm = T),rowMeans(NTG.v.hemi<0,na.rm=T)),1,min)

# compare G20 and NTG vulnerability values
NTG.v.hemi <- do.call(cbind,lapply(results.NTG,function(R) hemi.average(R$resid[,MPI.names[-1]])))
G20.v.hemi <- do.call(cbind,lapply(results.G20,function(R) hemi.average(R$resid[,MPI.names[-1]])))
NTG.minus.G20.v.hemi <- NTG.v.hemi - G20.v.hemi
p.v0.hemi <- 2*apply(cbind(rowMeans(NTG.minus.G20.v.hemi>0,na.rm = T),rowMeans(NTG.minus.G20.v.hemi<0,na.rm=T)),1,min)
p.lab <- setNames(rep('',n.regions.CNDR/2),names(p.v0.hemi))
p.lab[p.v0.hemi < 0.05] <- '*'
p.lab[p.v0.hemi < 0.01] <- '**'
p.lab[p.v0.hemi == 0] <- '***'
df.plt <- data.frame(diff=rowMeans(NTG.minus.G20.v.hemi),name=rownames(NTG.minus.G20.v.hemi),p=p.lab,stringsAsFactors = F)
df.plt <- df.plt[order(abs(df.plt$diff),decreasing = T)[1:30],]
p.v <- ggplot(data=df.plt) + geom_col(aes(x=name,y=diff),fill="#00A08A",alpha=0.5) +
  geom_text(aes(x=name,label=p,y=Inf),vjust=1,size=1.5)+
  scale_x_discrete(limits=df.plt$name[order(abs(df.plt$diff),decreasing = T)]) +
  ylab('NTG - G20 Vulnerability') + xlab('') + theme_classic() +
  theme(text=element_text(size=8),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
ggsave(p.v,filename = paste(savedir,'NTGvsG20BootstrapVulnerability_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 9,useDingbats=FALSE)
write.table(x='*, p<0.05; **, p < 0.01; ***, p < 0.002',file = paste0(savedir,'NTGvsG20BootstrapVulnerability_CNDRSpace_pkey.txt'))
