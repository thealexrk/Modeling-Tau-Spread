# plot connectivity from injection site vs pathology

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'connectivity_vs_path/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

# get mean pathology for each time point
tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
dimnames(W) <- list(region.names,region.names)

#########################################################################################
### for given injection site(s), extract connectivity matrix row/column in CNDR space ###
#########################################################################################

conn.inj.retro <- do.call(cbind,lapply(injection.site, function(s) map.ABA.to.CNDR(X.ABA = W[,s],CNDR.names = path.names,ABA.to.CNDR.key)))
colnames(conn.inj.retro) <- injection.site
conn.inj.antero <- do.call(cbind,lapply(injection.site, function(s) map.ABA.to.CNDR(X.ABA = W[s,],CNDR.names = path.names,ABA.to.CNDR.key)))
colnames(conn.inj.antero) <- injection.site

write.csv(conn.inj.retro,paste0(savedir,'RetrogradeConnectivityToInjectionSites.csv'))
write.csv(conn.inj.antero,paste0(savedir,'AnterogradeConnectivityToInjectionSites.csv'))

###############################################
### compare mean connectivity and pathology ###
###############################################

plot.col <- wes_palettes$Darjeeling1[4]

conn.retro.plot <- log(rowMeans(conn.inj.retro),base=10)
masks <- lapply(1:length(tps), function(t.) Grp.mean[[t.]] >0 & conn.retro.plot!=-Inf)
p.list <- lapply(1:length(tps), function(t.) p.xy(y = log(Grp.mean[[t.]][masks[[t.]]],base=10),x=conn.retro.plot[masks[[t.]]],xlab = 'log Retrograde Connectivity\nto Injection Site(s)',
                                                  ylab = paste0('log(',grp,' Pathology)'),ttl = paste(tps[t.],'MPI'),col = plot.col,alpha=0.8)+
                                                  theme(axis.title = element_text(size=6)))
p.all <- plot_grid(plotlist = p.list,nrow=1)
ggsave(p.all,filename = paste(savedir,grp,'_RetrogradeConnectivityToInjectionSites_vs_Pathology.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))

conn.antero.plot <- log(rowMeans(conn.inj.antero),base=10)
masks <- lapply(1:length(tps), function(t.) Grp.mean[[t.]] >0 & conn.antero.plot!=-Inf)
p.list <- lapply(1:length(tps), function(t.) p.xy(y = log(Grp.mean[[t.]][masks[[t.]]],base=10),x=conn.antero.plot[masks[[t.]]],xlab = 'log Anterograde Connectivity\nto Injection Site(s)',
                                                  ylab = paste0('log(',grp,') Pathology'),ttl = paste(tps[t.],'MPI'),col = plot.col,alpha=0.8)+
                                                  theme(axis.title = element_text(size=6)))
p.all <- plot_grid(plotlist = p.list,nrow=1)
ggsave(p.all,filename = paste(savedir,grp,'_AnterogradeConnectivityToInjectionSites_vs_Pathology.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps))
