# ask whether vulnerability values are consistent across time and hemispheres

#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe','injection.site')))
basedir <- params$basedir
setwd(basedir)
diffmodel.savedir <- paste(params$opdir,'diffmodel/bidirectional_onelm/',paste0(injection.site,collapse='-'),'/',sep='')
savedir <- paste(params$opdir,'diffmodel/vuln_time_hemi/',paste0(injection.site,collapse='-'),'/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
source('code/misc/plottingfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

##############################
### load path data for NTG ###
##############################

tps <- params$tps
Mice <- lapply(tps, function(tp) path.data[path.data$Condition == grp & path.data$Month == tp,path.names])
Grp.mean <- lapply(Mice,function(X) colMeans(X,na.rm = T))

####################################
### load vulnerability for group ###
####################################

vuln <- read.csv(file=paste0(diffmodel.savedir,grp,'vulnerability_bidirectional.csv'),check.names = F,row.names = 1)
vuln.hemi <- hemi.average(vuln)
colnames(vuln.hemi) <- 'Time-Hemi Average'

#################################################
### compare vulnerability over time and hemis ###
#################################################

vuln <- vuln[,which(colnames(vuln) != 'Average')] # remove average vulnerability column
list[vuln.i,vuln.c] <- hemi.split(vuln,rename=TRUE) # split vulnerability into two dataframes

colnames(vuln.i) <- paste(colnames(vuln.i),'ipsi')
colnames(vuln.c) <- paste(colnames(vuln.c),'contra')

if(identical(rownames(vuln.i),rownames(vuln.c)) & identical(rownames(vuln.i),rownames(vuln.hemi))){
  df.plt <- cbind(vuln.i,vuln.c,Average=vuln.hemi)  
}

# first compare distribution of absolute value of vulnerability
df.flat <- collapse.columns(df.plt,cnames=c(paste(tps,'MPI ipsi'),paste(tps,'MPI contra')))
df.flat$group <- substr(df.flat$names,1,5)
df.flat$names <- ifelse(grepl('ipsi',df.flat$names),yes='ipsi',no='contra')

pvals <- sapply(tps, function(t.) wilcox.test(x=df.plt[,paste(t.,'MPI ipsi')],y=df.plt[,paste(t.,'MPI contra')])$p.value)
plabs <- paste('p =',signif(pvals,2))
p <- ggplot(df.flat) + geom_jitter(aes(x=group,y=values,color=names),alpha=0.5,stroke=0,position=position_jitterdodge()) +
  annotate(geom='text',x=paste(tps,'MPI'),y=Inf,label=plabs,size=2.5,vjust=1)+ ylab('Vulnerability') + xlab('')+
  scale_color_brewer(palette='Set1') + theme_classic() + theme(text=element_text(size=8),legend.title = element_blank())
ggsave(p,filename = paste(savedir,grp,'CompareVulnerabilityBetweenHemispheres.pdf',sep=''),
       units = 'cm',height = 9,width = 10,useDingbats=FALSE)

# now plot spatial similarity of vulnerability
r.mat <- cor(df.plt,use = 'pairwise.complete.obs')
melt_mat <- melt(t(r.mat))
melt_mat$Var1 <- as.character(melt_mat$Var1)
melt_mat$value <- signif(melt_mat$value,2)
p <- imagesc(r.mat,cmap = 'redblue',clim = c(-1,1),caxis_name = 'r',ttl='Similarity of Model Residuals Across Time and Hemisphere') + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_text(data = melt_mat, aes(x=Var1, y=Var2, label=value),size=2.5) + nice_cbar('right') + coord_equal()

ggsave(p,filename = paste(savedir,grp,'VulnerabilityByTimeAndHemisphere.pdf',sep=''),
       units = 'cm',height = 9,width = 10)

## exclude 1 MPI
vuln.i <- vuln.i[,paste(c(3,6,9),'MPI ipsi')]
vuln.c <- vuln.c[,paste(c(3,6,9),'MPI contra')]
vuln.hemi <- hemi.average(vuln[,paste(c(3,6,9),'MPI')])
colnames(vuln.hemi) <- 'Time-Hemi Average'

if(identical(rownames(vuln.i),rownames(vuln.c)) & identical(rownames(vuln.i),rownames(vuln.hemi))){
  df.plt <- cbind(vuln.i,vuln.c,Average=vuln.hemi)  
}

r.mat <- cor(df.plt,use = 'pairwise.complete.obs')
melt_mat <- melt(t(r.mat))
melt_mat$Var1 <- as.character(melt_mat$Var1)
melt_mat$value <- signif(melt_mat$value,2)
p <- imagesc(r.mat,cmap = 'redblue',clim = c(-1,1),caxis_name = 'r',ttl='Similarity of Model Residuals Across Time and Hemisphere') + 
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_text(data = melt_mat, aes(x=Var1, y=Var2, label=value),size=2.5) + nice_cbar('right') + coord_equal()

ggsave(p,filename = paste(savedir,grp,'VulnerabilityByTimeAndHemisphere_Exclude1MPI.pdf',sep=''),
       units = 'cm',height = 9,width = 10)

