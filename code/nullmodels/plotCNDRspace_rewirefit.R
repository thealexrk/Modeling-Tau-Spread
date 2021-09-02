#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','injection.site')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'nullmodels/rewire/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
null.names <- c('InDegreePreserved','OutDegreePreserved')
col <- getGroupColors(params$grps)[grp] # get color for plots

for(null.name in null.names){
  
	load(paste(savedir,grp,'CNDRSpace',null.name,'BidirectionalOptim_data.RData',sep=''))
	tps <- params$tps

	# exclude regions with 0 pathology at each time point for purposes of computing fit
	lapply(1:length(tps), function(t) paste0(t,' MPI: ',sum(df[[t]]$path != -Inf & !is.na(df[[t]]$pred.retro)),'/',nrow(df[[t]]),' regions left'))
	
	# plot for each time point, using p.xy function 
	p <- lapply(1:length(tps), function(t) 
  p.xy(x=m[[t]]$fitted.values,y=m[[t]]$model$path,ylab='Actual',xlab='Predicted',
       ttl=paste0(grp,': ',tps[t],' MPI'),col=col,alpha=0.7))
	p <- plot_grid(plotlist=p,align='hv',nrow=1)

	ggsave(p,filename = paste(savedir,grp,'CNDRSpaceFit_',null.name,'.pdf',sep=''),
	       units = 'cm',height = 3.75,width = 3.75*length(tps),useDingbats=FALSE)
}