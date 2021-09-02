#################
### Load data ###
#################
grp <- 'IgG1'
rm(list=setdiff(ls(),c('params','grp')))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'diffmodel/',sep='')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
source('code/misc/miscfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
#load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression

# get mean pathology for only time point, 1 month
tp <- 1
Mice <- path.data[path.data$Condition == grp,-1]
Grp.mean <- colMeans(Mice,na.rm = T)

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W

ABA.missing <- data.frame(region_names=names(CNDR.to.ABA.key)[sapply(CNDR.to.ABA.key,length)==0])
write.csv(ABA.missing,file = paste0(params$opdir,'processed/ABARegionsWithNoCNDRMatch.csv'))

################################
### Test: ABA -> CNDR -> ABA ###
################################

nperms <- 1000
# save: distribution of r values for retrieving random gaussians (1), excluding identical regions (2), number of identical regions (3)
results <- data.frame(r.dist=rep(NA,nperms),r.dist.nonidentical=rep(NA,nperms),num.identical=rep(NA,nperms))

for(p in 1:nperms){
  print(paste('perm',p))
  X.ABA.groundtruth <- rnorm(n.regions.ABA) # define ground truth
  names(X.ABA.groundtruth) <- region.names
  CNDR.names <- path.names
  X.CNDR <- quiet(map.ABA.to.CNDR(X.ABA.groundtruth,CNDR.names,ABA.to.CNDR.key)) # suppress output for now
  
  ABA.names <- region.names
  X.ABA.retrieve <- quiet(map.CNDR.to.ABA(X.CNDR,ABA.names,CNDR.to.ABA.key)) # suppress output for now
  results$r.dist[p] <- cor(X.ABA.groundtruth,X.ABA.retrieve,use='pairwise.complete.obs')
  non.identical.idx <- which(!(X.ABA.groundtruth==X.ABA.retrieve))
  results$num.identical <- sum(X.ABA.groundtruth==X.ABA.retrieve,na.rm=T)
  results$r.dist.nonidentical[p] <- cor(X.ABA.groundtruth[non.identical.idx],X.ABA.retrieve[non.identical.idx],use='pairwise.complete.obs')
}

tsz <- 2
c1 <- wes_palettes$Rushmore1[2]
c2 <- wes_palettes$Rushmore1[3]
# plot a single test for visualization
df <- data.frame(x= X.ABA.groundtruth[non.identical.idx], y = X.ABA.retrieve[non.identical.idx])
n.missing <- sum(is.na(X.ABA.retrieve))
n.identical <- unique(results$num.identical)
r <- signif(cor(df$x,df$y,use = 'pairwise.complete.obs'),2)
p1 <- ggplot(data=df) + geom_point(aes(x=x,y=y),alpha= 0.8,stroke=0,color = c1) +
  xlab('Ground Truth') + ylab('Non-identically Retrieved') + ggtitle('ABA to CNDR to ABA') +
  annotate(geom = 'text',label=paste0(n.missing,'/',n.regions.ABA,' missing \n',n.identical,' identical'),x=max(df$x),y=min(df$y,na.rm = T),hjust=1,vjust=0,size=tsz) +
  annotate(geom = 'text',label=paste0('r = ',r),x=min(df$x),y=max(df$y,na.rm = T),hjust=0,size=tsz) +
  theme_classic() + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))

# plot distribution of retrieved correlation values overall and 
df <- data.frame(x=results$r.dist.nonidentical,x2=results$r.dist)
p2 <- ggplot(data=df) + geom_histogram(aes(x=x),alpha = 0.3,fill = c2) + xlab('r(Ground Truth, Non-identically Retrieved)') +
  scale_x_continuous(limits=c(0,1)) + ylab('Count') +  ggtitle('ABA to CNDR to ABA') + 
  theme_classic() + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
p2
# combine all
p <- plot_grid(plotlist = list(p1,p2),align = 'hv',nrow=1)
ggsave(p,filename = paste0(params$opdir,'processed/ABAtoCNDRtoABA.pdf'),units = 'in',height = 2,width = 6)

#################################
### Test: CNDR -> ABA -> CNDR ###
#################################

# save: distribution of r values for retrieving random gaussians (1), excluding identical regions (2), number of identical regions (3)
results <- data.frame(r.dist=rep(NA,nperms),r.dist.nonidentical=rep(NA,nperms),num.identical=rep(NA,nperms))

for(p in 1:nperms){
  print(paste('perm',p))
  X.CNDR.groundtruth <- rnorm(n.regions.CNDR) # define ground truth
  names(X.CNDR.groundtruth) <- path.names
  ABA.names <- region.names
  X.ABA <- quiet(map.CNDR.to.ABA(X.CNDR.groundtruth,ABA.names,CNDR.to.ABA.key)) # suppress output for now
  
  CNDR.names <- path.names
  X.CNDR.retrieve <- quiet(map.ABA.to.CNDR(X.ABA,CNDR.names,ABA.to.CNDR.key)) # suppress output for now
  results$r.dist[p] <- cor(X.CNDR.groundtruth,X.CNDR.retrieve,use='pairwise.complete.obs')
  non.identical.idx <- which(!(X.CNDR.groundtruth==X.CNDR.retrieve))
  results$num.identical <- sum(X.CNDR.groundtruth==X.CNDR.retrieve,na.rm=T)
  results$r.dist.nonidentical[p] <- cor(X.CNDR.groundtruth[non.identical.idx],X.CNDR.retrieve[non.identical.idx],use='pairwise.complete.obs')
}

c1 <- wes_palettes$Rushmore1[2]
c2 <- wes_palettes$Rushmore1[3]
tsz <- 2
# plot a single test for visualization
df <- data.frame(x= X.CNDR.groundtruth[non.identical.idx], y = X.CNDR.retrieve[non.identical.idx])
n.missing <- sum(is.na(X.CNDR.retrieve))
n.identical <- unique(results$num.identical)
r <- signif(cor(df$x,df$y,use = 'pairwise.complete.obs'),2)
p1 <- ggplot(data=df) + geom_point(aes(x=x,y=y),alpha= 0.8,stroke=0,color = c1) +
  xlab('Ground Truth') + ylab('Non-identically Retrieved') + ggtitle('CNDR to ABA to CNDR') +
  annotate(geom = 'text',label=paste0(n.missing,'/',n.regions.CNDR,' missing \n',n.identical,' identical'),x=max(df$x),y=min(df$y,na.rm = T),hjust=1,vjust=0,size=tsz) +
  annotate(geom = 'text',label=paste0('r = ',r),x=min(df$x),y=max(df$y,na.rm = T),hjust=0,size=tsz) +
  theme_classic() + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))

# plot distribution of retrieved correlation values overall and 
df <- data.frame(x=results$r.dist.nonidentical,x2=results$r.dist)
p2 <- ggplot(data=df) + geom_histogram(aes(x=x),alpha = 0.3,fill = c2) + xlab('r(Ground Truth, Non-identically Retrieved)') +
  scale_x_continuous(limits=c(0,1)) + ylab('Count') +  ggtitle('CNDR to ABA to CNDR') + 
  theme_classic() + theme(text=element_text(size=8),plot.title = element_text(hjust=0.5))
p2
# combine all
p <- plot_grid(plotlist = list(p1,p2),align = 'hv',nrow=1)
ggsave(p,filename = paste0(params$opdir,'processed/CNDRtoABAtoCNDR.pdf'),units = 'in',height = 2,width = 6)
