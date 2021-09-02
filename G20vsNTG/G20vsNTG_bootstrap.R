#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','injection.site')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/bidirectional_bootstrap/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/optimfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

###################################################
### load bootstrap data for bidirectional model ###
###################################################

load(file = paste(savedir,'NTGCNDRSpaceBidirectionalOptim_Bootstrap.RData',sep=''))
results.NTG <- results

load(file = paste(savedir,'G20CNDRSpaceBidirectionalOptim_Bootstrap.RData',sep=''))
results.G20 <- results
rm(results)

MPI.names <- paste(params$tps,'MPI')
group.colors <- getGroupColors(params$grps)

# 1. compare model fits between g20 and NTG
cfg <- list(train.fits.mse='MSE',train.fits.r='Pearson r') # loop through pearson r and MSE
p.fit <- list()
for(f.met in names(cfg)){
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
  subset.list <- lapply(c(Anterograde='Anterograde',Retrograde='Retrograde'), function(n) df.plt$values[df.plt$names == n && df.plt$group==grp])
  subset.stats <- lapply(subset.list, function(X) c(mean=mean(X),quantile(X,c(0.025,.975)),median=median(X),cv=sd(X)/mean(X)))
  sink(paste0(savedir,grp,"TimeConstantStats.txt"))
  print(subset.stats)
  sink()
}
p.lab <- sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0,df.plt$values[df.plt$names==d & df.plt$group == 'NTG']-df.plt$values[df.plt$names==d & df.plt$group == 'G20']))
p.lab <- paste0('p = ',signif(p.lab,2))
df.plt$group <- factor(df.plt$group,levels = c(params$grps),ordered=T)
p.c <- ggplot(df.plt) + geom_boxplot(aes(x=names,y=values,fill=group),size=0.25,outlier.size = 0.25) + theme_classic() +
  annotate(geom='text',x=c('Anterograde','Retrograde'),y=Inf,label=p.lab,vjust=1,size=2)+
  ylab('Diffusivity Constant') + xlab('') + scale_y_continuous(limits=c(0,NA))+scale_fill_manual(limits=params$grps,values =group.colors,name='')+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
p.c
ggsave(p.c,filename = paste(savedir,'NTGvsG20TimeConstants_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 4,width = 5,useDingbats=FALSE)

# 3. compare anterograde and retrograde betas between groups
cfg <- list(NTG=results.NTG,G20=results.G20) # loop through group results
for(grp in params$grps){
  ant.coefs <-as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) setNames(sapply(R$m.train, function(m.i) coefficients(lm.beta(m.i))['pred.antero']),MPI.names))))
  ret.coefs <-as.data.frame(do.call(rbind,lapply(cfg[[grp]], function(R) setNames(sapply(R$m.train, function(m.i) coefficients(lm.beta(m.i))['pred.retro']),MPI.names))))
  df.plt <- rbind(cbind(ant.coefs,OI='Anterograde',stringsAsFactors=F),cbind(ret.coefs,OI='Retrograde',stringsAsFactors=F))
  colnames(df.plt) <- c(MPI.names,'OI')
  df.plt <- collapse.columns(df.plt,cnames = MPI.names,groupby = 'OI')
  df.plt$Group <- grp
  p.avr <- sapply(MPI.names, function(M) pval.2tail.np(0,df.plt$values[df.plt$names==M & df.plt$group=='Anterograde'] - df.plt$values[df.plt$names==M & df.plt$group=='Retrograde']))
  p.avr <- p.avr
  p.lab.avr <- sapply(p.avr, function(P) paste0('p = ',signif(P,2)))
  p.lab.avr[p.avr ==0] <- paste0('p < ',signif(1/length(cfg[[grp]]),2)) # don't say p = 0, say < 1/nboots
  # prep data for making line plots with error bars -- i hate this so much and i don't think i'll use it but leaving here for now
  line.row.names <- as.vector(sapply(MPI.names, function(m) sapply(c('Anterograde','Retrograde'),function(d) paste0(m,d))))
  df.line <- matrix(NA, nrow=length(line.row.names),ncol=5,dimnames=list(line.row.names,c('Month','Direction','Mean','Upper','Lower')))
  for(tp in MPI.names){
    for(d in c('Anterograde','Retrograde')){
      rn <- paste0(tp,d)
      df.line[rn,'Direction'] <- d
      df.line[rn,'Month'] <- as.numeric(substr(tp,1,1))
      df.line[rn,'Mean'] <- mean(df.plt$values[df.plt$names==tp & df.plt$group==d])
      df.line[rn,'Upper'] <- quantile(df.plt$values[df.plt$names==tp & df.plt$group==d],0.975)
      df.line[rn,'Lower'] <- quantile(df.plt$values[df.plt$names==tp & df.plt$group==d],0.025)
    }
  }
  df.line <- as.data.frame(df.line,stringsAsFactors=F)
  for(col in c('Mean','Upper','Lower')){df.line[,col] <- as.numeric(df.line[,col])}
  df.line$Group <- grp
  cfg[[grp]] <- list(df=df.plt,p.avr=p.lab.avr,grp=grp,df.line=df.line) # store plot ready df in looped list
}

df.plt <- do.call(rbind,lapply(cfg,function(X) X$df)) # vertically concatenate df across G20 and NTG then plot
# p-value comparing ant vs. ret for each group
p.avr.save <- do.call(rbind,lapply(cfg, function(X) c(X$p.avr,Group=X$grp,Pcrit=0.05/length(X$p.avr))))
write.csv(x=p.avr.save,file=paste0(savedir,'NTGandG20_AnteroVsRetroBetas_Stats.csv'),row.names = F)
# compare the two groups to each other
p.G20.v.NTG <- sapply(MPI.names, function(M)  # loop through months
  # isolate ant and ret for each month and compare between groups
                sapply(c('Anterograde','Retrograde'), function(d) pval.2tail.np(0, # loop through direction
                  cfg[['NTG']]$df$values[cfg[['NTG']]$df$names==M & cfg[['NTG']]$df$group==d] - # compare group diff to 0
                    cfg[['G20']]$df$values[cfg[['G20']]$df$names==M & cfg[['G20']]$df$group==d])))
p.G20.v.NTG <- list.posthoc.correct(p.G20.v.NTG,'bonferroni') # bonferroni correct
p.G20.v.NTG.lab <- as.data.frame(sapply(as.data.frame(p.G20.v.NTG),function(x) paste('p =',signif(x,2))))
p.G20.v.NTG.lab$group <- c('Anterograde','Retrograde')
df.p <- collapse.columns(p.G20.v.NTG.lab,cnames=MPI.names,groupby='group')
df.plt$Group <- factor(df.plt$Group,levels = params$grps,ordered=T)
p <- ggplot() + geom_boxplot(data=df.plt,aes(x=names,y=values,fill=Group),size=0.25,outlier.size=0.25) + 
  theme_classic() + facet_wrap(~group)+
  geom_text(data=df.p,aes(x=names,label=values,y=Inf),vjust=1,size=1.8)+
  ylab('Standardized Beta') + xlab('') +scale_fill_manual(limits=params$grps,values =group.colors,name='')+
  theme(text=element_text(size=8),legend.key.size = unit(0.1,'cm'),legend.box.margin = ggplot2::margin(t = 0, unit='cm'),axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(p,filename = paste(savedir,'NTGvsG20AnterogradeRetrogradeBetas_Boxplot_CNDRSpace.pdf',sep=''),
       units = 'cm',height = 6,width = 9,useDingbats=FALSE)

# 4. bootstrap NTG vulnerability values

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
