#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/euclidean/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(savedir,grp,'CNDRSpaceEuclideanDistanceFit_data.RData',sep=''))
tps <- params$tps
col <- getGroupColors(params$grps)[grp] # get color for plots

# exclude regions with 0 pathology at each time point for purposes of computing fit
lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred)),'/',nrow(df[[t]]),' regions left'))
# plot for each time point, using p.xy function 
p <- lapply(1:length(tps), function(t) 
  p.xy(x=df[[t]]$pred,y=df[[t]]$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col=col,alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_Euclidean.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)

# plot fit and label each region
geom_label_repel(aes(label = Region),
                 point.padding = 0.5,
                 segment.color = 'grey50')
p <- list()
for(t. in 1:length(tps)){
  t.path.names <- substr(rownames(df[[t.]]),1,nchar(rownames(df[[t.]])))
  df.tmp <- data.frame(x=df[[t.]]$pred,y=df[[t.]]$path,lab=t.path.names)
  df.tmp$z <- scale(df.tmp$x,center=T)
  #df.tmp <- df.tmp[substr(rownames(df.tmp),1,1) == 'i',]
  df.tmp <- df.tmp[abs(df.tmp$z) >1.5,]
  p[[t.]] <- p.xy(x=df[[t.]]$pred,y=df[[t.]]$path,ylab='Actual',xlab='Predicted',
                  ttl=paste0(grp,': ',tps[t.],' MPI'),col=col,alpha=0.7) + geom_label_repel(data=df.tmp,aes(x=x,y=y,label=lab),
                                                                                            point.padding = 0.5,label.size=0.05,size=1,label.r=0.1,
                                                                                            segment.color = 'grey50')
}

p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_Euclidean_label_repel.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)

# same plot as above excluding injection sites. i don't know why these are outliers consistently.
# it could be because the distance matrices are fully dense so pathology always gets back there?

injection.site.CNDR <- unname(params$injection.site.CNDR[injection.site]) # convert injection site from ABA to CNDR
df.excl <- lapply(df, function(X) X[path.names[!path.names %in% injection.site.CNDR],])

p <- lapply(1:length(tps), function(t) 
  p.xy(x=df.excl[[t]]$pred,y=df.excl[[t]]$path,ylab='Log(Pathology)',xlab='Log(Predicted)',
       ttl=paste0(grp,': ',tps[t],' MPI'),col=col,alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)

ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_Euclidean_ExcludeInjectionSites.pdf',sep=''),
       units = 'cm',height = 4.35,width = 4.35*length(tps),useDingbats=FALSE)

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
       ttl=paste0(grp,': ',tps[t],' MPI'),col=col,alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceVulnerability_Euclidean_vs',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)

##############################
### Include genes in model ###
##############################

m <- lapply(df.gene, function(df.i) m <- lm(path~.,data=inf.nan.mask(df.i)))
vulnerability <- lapply(m,residuals)
pred <- lapply(m, function(m.i) m.i$fitted.values)

p <- lapply(1:length(tps), function(t) 
  p.xy(x=m[[t]]$fitted.values,y=m[[t]]$model$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col=col,alpha=0.7))
p <- plot_grid(plotlist=p,align='hv',nrow=1)
ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_Euclidean+',goi,probe,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)

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