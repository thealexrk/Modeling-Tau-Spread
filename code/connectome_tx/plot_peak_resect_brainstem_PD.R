# how does deleting connections affect overall pathology and model fit

#################
### Load data ###
#################

grp <- 'NTG'

rm(list=setdiff(ls(),c('params','grp')))
print(grp)
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'connectome_tx/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')

ant.ret <- 'retro'
load(file = paste(savedir,grp,'FiberResectionPeakAUC',ant.ret,'.RData',sep=''))
tps <- params$tps
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

W <- readMat(paste0(params$opdir,'processed/W.mat'))$W
in.strength <- colSums(W)
names(in.strength) <- region.names
out.strength <- rowSums(W)
names(out.strength) <- region.names

############
### Plot ###
############

# ask whether for multiple brain stem sites, there are regions that have a large effect on impulse response of
# pathology curves, out of proportion to what you would expect based on injection site.

roi.lab.static <- c('iSNc','iSNr') # label the iSN in all plots
ptcol <- "#7294D4"
p.aupc.list <- p.aupc.list.no.out <- p.aupc.instrength.list.no.out<- p.peak.list <- p.aupc.hist.list <- list()
for(injection.site.name in names(brainstem.PD.sites)[-5]){
  peak.diff.matrix <- deleted.sumstats[[injection.site.name]]$peak.diff.matrix
  conn.to.inj <- deleted.sumstats[[injection.site.name]]$conn.to.inj
  if(length(dim(conn.to.inj))>1){conn.to.inj <- rowMeans(conn.to.inj)}
  names(conn.to.inj) <- region.names[which(!region.names %in% brainstem.PD.sites[[injection.site.name]])]
  aupc.del <- deleted.sumstats[[injection.site.name]]$aupc.del
  # remove any regions with 0 connectivity to injection site b/c deletion would have no effect
  remove.no.conn <- conn.to.inj > 0
  conn.to.inj <- conn.to.inj[remove.no.conn]
  aupc.del <- aupc.del[remove.no.conn] - aupc.baseline[[injection.site.name]]
  in.strength.inj <- in.strength[names(aupc.del)]
  total.delta.peak <- rowSums(peak.diff.matrix[remove.no.conn,]) # sum up the change in peak times over all nodes, for each deleted connection
  
  # roi.lab <- c(roi.lab.static,names(total.delta.peak)[which.min(total.delta.peak)]) # label one region in all (iSN) and then the most affected region
  # p.peak.list[[injection.site.name]] <- p.xy(conn.to.inj,total.delta.peak,xlab = 'Connection Strength\nwith Injection Site',ylab='Total Change\nin Peak Times',
  #           ttl = paste(injection.site.name,'Seed'),col=ptcol) +
  #   annotate(geom='text',x=conn.to.inj[roi.lab],y=total.delta.peak[roi.lab],label=roi.lab,hjust=1.2,color='red',alpha=0.7,size=2)
  # ggsave(p.peak.list[[injection.site]],filename = paste(savedir,grp,'FiberResectionTotalPeakTimeChangevsConnectionStrengthToSeed',injection.site,'.pdf',sep=''),
  #        units = 'cm',height = 4,width = 4)
  
  roi.lab <- c(roi.lab.static,names(aupc.del)[which.min(aupc.del)]) # label one region in all (iSN) and then the most affected region
  p.aupc.list[[injection.site.name]] <- p.xy(conn.to.inj,aupc.del,xlab = 'Connectivity with Injection Site\n',ylab=expression(Delta*" AUPC"),
                                        ttl = paste(injection.site.name,'Seed'),col=ptcol) #+
    #annotate(geom='text',x=conn.to.inj[roi.lab],y=aupc.del[roi.lab],label=roi.lab,hjust=1.2,color='red',alpha=0.7,size=2)
  p.aupc.list.no.out[[injection.site.name]] <- p.xy(conn.to.inj[!outlier.mask(aupc.del) & !outlier.mask(conn.to.inj)],aupc.del[!outlier.mask(aupc.del)& !outlier.mask(conn.to.inj)],
                                                    xlab = 'Connectivity with Injection Site\n',ylab=expression(Delta*" AUPC"),
                                             ttl = paste(injection.site.name,'Seed'),col=ptcol)
  p.aupc.instrength.list.no.out[[injection.site.name]] <- p.xy(in.strength.inj,aupc.del,
                                                    xlab = 'In-Strength of Deleted Edge\n',ylab=expression(Delta*" AUPC"),
                                                    ttl = paste(injection.site.name,'Seed'),col=ptcol)
 p.aupc.hist.list[[injection.site.name]] <- local({ 
   aupc.del <- aupc.del; injection.site.name <- injection.site.name
   aupc.out <- sort(aupc.del[outlier.mask(aupc.del)]) # get outliers in area under the path curve
   n.out <- min(length(aupc.out),7); # only show top n outliers if existing
   xmid <- min(aupc.del) + 0.5*(max(aupc.del)-min(aupc.del))
   col.scl <- rev(colorRampPalette(brewer.pal(9, 'Reds'))(n.out+3))
 p <- ggplot() + geom_histogram(aes(x=aupc.del),fill=ptcol) + 
   geom_histogram(aes(x=aupc.out),fill='red')
  for(o in 1:n.out){p <- p + annotate(geom='text',x=xmid,y=Inf,vjust=1.2*o,label=names(aupc.out[o]),color=col.scl[o],size=2.5)}
  p <- p + theme_classic() +
   xlab(expression(Delta*" AUPC"))+ggtitle(paste(injection.site.name,'Seed')) + theme(plot.title = element_text(hjust=0.5,face='bold'),text=element_text(size=8))
  print(p)
  })
}

p.all <- plot_grid(plotlist = p.aupc.list,nrow = 2,ncol= 3)
ggsave(p.all,filename = paste(savedir,grp,'FiberResectionAUPCvsConnectionStrengthToSeed_AllSeeds.pdf',sep=''),
       units = 'cm',height = 12,width = 18)
p.all <- plot_grid(plotlist = p.aupc.list.no.out,nrow = 2,ncol= 3)
ggsave(p.all,filename = paste(savedir,grp,'FiberResectionAUPCvsConnectionStrengthToSeed_AllSeedsNoOutliers.pdf',sep=''),
       units = 'cm',height = 12,width = 18)
p.all <- plot_grid(plotlist = p.aupc.instrength.list.no.out,nrow = 2,ncol= 3)
ggsave(p.all,filename = paste(savedir,grp,'FiberResectionAUPCvsDeletedEdgeInStrength_AllSeedsNoOutliers.pdf',sep=''),
       units = 'cm',height = 12,width = 18)
# p.all <- plot_grid(plotlist = p.peak.list,nrow = 2,ncol=3)
# ggsave(p.all,filename = paste(savedir,grp,'FiberResectionTotalPeakTimeChangevsConnectionStrengthToSeed_AllSeeds.pdf',sep=''),
#        units = 'cm',height = 12,width = 18)

p.all <- plot_grid(plotlist = p.aupc.hist.list,nrow = 2,ncol=3)
ggsave(p.all,filename = paste(savedir,grp,'FiberResectionAUPCOutlierHistogram_AllSeeds.pdf',sep=''),
       units = 'cm',height = 12,width = 18)

zscore <- function(x) scale(x,center =T)
df.list <- lapply(deleted.sumstats, function(X) data.frame(AUPC=zscore(X$aupc.del[X$conn.to.inj>0]),Conn=zscore(X$conn.to.inj[X$conn.to.inj>0])))
df.all <- do.call(rbind,df.list)
df.all <- df.all[!outlier.mask(df.all$Conn) & !outlier.mask(df.all$AUPC),]
# p.xy(x = df.all$Conn,df.all$AUPC,xlab='Connectivity to Injection',ylab=expression(Delta*"AUPC"))
# df.list <- lapply(names(df.list), function(inj) cbind(df.list,Inj=))