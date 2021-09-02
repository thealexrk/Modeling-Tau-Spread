# print which regions match and don't match between ABA and CNDR mappings

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

X.ABA <- c(1:length(region.names))
names(X.ABA) <- region.names
CNDR.names <- path.names
X.CNDR <- map.ABA.to.CNDR(X.ABA,CNDR.names,ABA.to.CNDR.key)

X.CNDR <- c(1:length(path.names))
names(X.CNDR) <- path.names
ABA.names <- region.names
X.ABA <- map.CNDR.to.ABA(X.CNDR,ABA.names,CNDR.to.ABA.key)

names(CNDR.to.ABA.key)[sapply(CNDR.to.ABA.key,length)==0]
