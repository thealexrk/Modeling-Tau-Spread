#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/retrograde/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(savedir,grp,'CNDRSpaceFit_data.RData',sep=''))
tps <- params$tps

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred,y=df[[t]]$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_basemodel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

c.rng <- seq(params$c.min,params$c.max,length.out = params$c.n)
p <- plot.Xt(Xt = t(Xt.sweep),t. = c.rng) + geom_vline(xintercept = c.Grp,color='grey50',linetype='dashed') +
  annotate(geom='text',x=c.Grp,y=max(Xt.sweep),hjust=-0.25,label='Optimal',size=2,color='grey50')+
  scale_color_manual(values = as.character(1:length(tps)),labels=paste(tps,'MPI'),name='') + 
  xlab('c') + ylab('Pearson r with\nPathology') +ggtitle(paste0(grp,': ',paste0(injection.site,collapse = '-'))) +
  theme_bw() + theme(text=element_text(size=8),plot.title = element_text(size=8,hjust=0.5),legend.key.size = unit(0.1,'cm'))
p
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceCSweepByTimePoint_basemodel.pdf',sep=''),
       units = 'cm',height = 4,width = 9)

###############################################
### Save predicted values and vulnerability ###
###############################################

m <- lapply(df, function(df.i) m <- lm(path~pred,data=inf.nan.mask(df.i))) # compute residuals
vulnerability <- lapply(m,residuals)
# initialize empty data frame with all CNDR regions and time points
df.vuln <- data.frame(matrix(ncol=length(tps),nrow=n.regions.CNDR, dimnames=list(path.names, paste(tps,'MPI'))),check.names = FALSE)
for(t in 1:length(tps)){
  regions.mpi <- names(vulnerability[[t]]) # get the regions that had non-zero pathology at each tp (so log can be computed)
  df.vuln[regions.mpi,paste(tps[t],'MPI')] <- vulnerability[[t]][regions.mpi]
}

write.csv(hemi.average(df.vuln),paste0(savedir,grp,'vulnerability_hemiaverage.csv'))
df.vuln <- cbind(df.vuln,Average=rowMeans(df.vuln,na.rm = T)) # compute average vulnerability across time points, ignoring regions with 0 path
write.csv(df.vuln,paste0(savedir,grp,'vulnerability.csv'))

df.pred <- do.call(cbind,lapply(df,function(X) X[,'pred',drop=FALSE])) # get log10 predicted pathology
write.csv(df.pred,paste0(savedir,grp,'log10predictedpath.csv'))

#####################
### Add gene data ###
#####################

# load gene data specified by goi, probe input variables
gene.exp <- read.csv(paste0(basedir,'data/aba/expression/',goi,probe,'ExpressionCNDR.csv'),row.names = 1)
df.gene <- lapply(df, function(x) merge(x,gene.exp,by=0)) # merge with pathology data frame
for(j in 1:length(df.gene)){rownames(df.gene[[j]]) <- df.gene[[j]]$Row.names} # remove row names
for(j in 1:length(df.gene)){df.gene[[j]] <- df.gene[[j]][,-1]}

##########################################################################
### show how gene expression correlates with residuals of spread model ###
##########################################################################

vuln.gene <- lapply(vulnerability, function(v) merge(as.data.frame(v),gene.exp,by=0))
p <- lapply(1:length(tps), function(t) 
  p.xy(x=vuln.gene[[t]][,paste0(goi,'_',probe)],y=vuln.gene[[t]]$v,ylab='Vulnerability',xlab=paste0(goi,' Expression'),
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceVulnerability_basemodel_vs',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

##############################
### Include genes in model ###
##############################

m <- lapply(df.gene, function(df.i) m <- lm(path~.,data=inf.nan.mask(df.i)))
vulnerability <- lapply(m,residuals)
pred <- lapply(m, function(m.i) m.i$fitted.values)

p <- lapply(1:length(tps), function(t) 
  p.xy(x=m[[t]]$fitted.values,y=m[[t]]$model$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col='#007257',alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_basemodel+',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

# concatenate predicted path (fitted vals, y-hat) and vulnerability (residuals, y-yhat) into one dataframe
# initialize empty data frame with all CNDR regions and time points
df.vuln <- df.pred <- data.frame(matrix(ncol=length(tps),nrow=n.regions.CNDR, dimnames=list(path.names, paste(tps,'MPI'))),check.names = FALSE)
for(t in 1:length(tps)){
  regions.mpi <- names(vulnerability[[t]]) # get the regions that had non-zero pathology at each tp (so log can be computed)
  df.vuln[regions.mpi,paste(tps[t],'MPI')] <- vulnerability[[t]][regions.mpi]
  df.pred[regions.mpi,paste(tps[t],'MPI')] <- pred[[t]][regions.mpi]
}

write.csv(hemi.average(df.vuln),paste0(savedir,grp,'vulnerability',goi,probe,'_hemiaverage.csv'))
df.vuln <- cbind(df.vuln,Average=rowMeans(df.vuln,na.rm = T)) # compute average vulnerability across time points, ignoring regions with 0 path
write.csv(df.vuln,paste0(savedir,grp,'vulnerability',goi,probe,'.csv'))
write.csv(df.pred,paste0(savedir,grp,'predictedpath',goi,probe,'.csv'))