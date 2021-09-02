#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','injection.site')))
basedir <- params$basedir
setwd(basedir)
diffmodel.savedir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
savedir <- paste(params$opdir,'G20vsNTG/path_vs_vuln/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

#####################################
### load path data for each group ###
#####################################

tps <- params$tps
Mice.NTG <- lapply(tps, function(tp) path.data[path.data$Condition == 'NTG' & path.data$Month == tp,path.names])
Grp.mean.NTG <- lapply(Mice.NTG,function(X) colMeans(X,na.rm = T))
Mice.G20 <- lapply(tps, function(tp) path.data[path.data$Condition == 'G20' & path.data$Month == tp,path.names])
Grp.mean.G20 <- lapply(Mice.G20,function(X) colMeans(X,na.rm = T))

##################################
### load vulnerability for NTG ###
##################################
NTG.vuln <- read.csv(file=paste0(diffmodel.savedir,'NTGvulnerability_bidirectional.csv'),check.names = F)
tp.excl <- '1 MPI' # exclude 1 MPI - see diffmodel/vuln_hemi_time.R
NTG.vuln.hemi <- t(read.csv(file=paste0(diffmodel.savedir,'NTGvulnerability_bidirectional_hemiaverage_exclude',tp.excl,'.csv'),row.names = 1))[1,]

#####################################################################################
### plot log(G20/NTG) vs. NTG vulnerability as in Fig 7d of Henderson et al. 2019 ###
#####################################################################################

#colorvis(wes_palettes$Zissou1)
plot.col <- wes_palettes$Zissou[2]

p.list <- lapply(1:length(tps), function(t.) p.xy(y = log(Grp.mean.G20[[t.]]/Grp.mean.NTG[[t.]],base=10),x=NTG.vuln[,paste(tps[t.],'MPI')],xlab = 'NTG Vulnerability',
                                       ylab = 'log(G2019S/NTG) Pathology',ttl = paste(tps[t.],'MPI'),col = plot.col,alpha=0.8))
p.all <- plot_grid(plotlist = p.list,nrow=1)
ggsave(p.all,filename = paste(savedir,'G2019-NTGvsNTGVulnerability_Bidirectional_TimeDependent.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)

######################################################################################
### plot log(G20/NTG) vs. NTG vulnerability with hemisphere averaged vulnerability ###
######################################################################################

G20.NTG.hemiaverage <- lapply(1:length(tps), function(t.) log(hemi.average(Grp.mean.G20[[t.]],TRUE)/hemi.average(Grp.mean.NTG[[t.]],TRUE),base=10))
p.list <- lapply(1:length(tps), function(t.) p.xy(y = G20.NTG.hemiaverage[[t.]],x=NTG.vuln.hemi,xlab = 'NTG Vulnerability',
                                                  ylab = 'log(G2019S/NTG) Pathology',ttl = paste(tps[t.],'MPI'),col = plot.col,alpha=0.8))
p.all <- plot_grid(plotlist = p.list,nrow=1)
ggsave(p.all,filename = paste(savedir,'G2019-NTGvsNTGVulnerability_Bidirectional_HemiTimeAverage_exclude',tp.excl,'.pdf',sep=''),
       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)
