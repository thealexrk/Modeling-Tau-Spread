#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)
source('code/misc/miscfxns.R')

path.data <- read.csv('data/PathData.csv',header = TRUE, check.names = FALSE)
connectivity.ipsi <- read.csv('data/Connectome_Ipsi.csv',row.names = 1,header = TRUE,check.names = FALSE)
connectivity.contra <- read.csv('data/Connectome_Contra.csv',row.names = 1,header = TRUE, check.names = FALSE)
ABA.to.CNDR <- read.csv('data/ABAtoCNDR.csv',header = TRUE, stringsAsFactors = F) # lists CNDR regions for each ABA region
CNDR.to.ABA <- read.csv('data/CNDRtoABA.csv',header = TRUE, stringsAsFactors = F) # lists ABA regions for each CNDR region

###################################
### Process connectivity matrix ###
###################################

conn.names.ipsi <- colnames(connectivity.ipsi)
conn.names.contra <- colnames(connectivity.contra)

# checks 
if(identical(colnames(connectivity.contra),rownames(connectivity.contra))){
  print('contra connectivity colnames and rownames equal')}
if(identical(colnames(connectivity.ipsi),rownames(connectivity.ipsi))){
  print('ipsi connectivity colnames and rownames equal')}

W <- rbind(cbind(connectivity.ipsi,connectivity.contra),cbind(connectivity.contra,connectivity.ipsi))
rownames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral regions
colnames(W) <- c(paste('i',rownames(connectivity.ipsi),sep=''), # add i to ipsilateral regions
                 paste('c',rownames(connectivity.contra),sep='')) # add c to contralateral regions

n.regions.ABA <- nrow(W)
n.regions.ABA.hemi <- n.regions.ABA/2
# check if connectivity was tiled into alternating blocks correctly

unit.test(all(W[(n.regions.ABA.hemi+1):n.regions.ABA,(n.regions.ABA.hemi+1):n.regions.ABA] == W[1:n.regions.ABA.hemi,1:n.regions.ABA.hemi]),
          'tiling on-diagonal blocks worked','tiling on-diagonal blocks failed') 
unit.test(all(W[1:n.regions.ABA.hemi,(n.regions.ABA.hemi+1):n.regions.ABA] == W[(n.regions.ABA.hemi+1):n.regions.ABA,1:n.regions.ABA.hemi]),
          'tiling off-diagonal blocks worked','tiling off-diagonal blocks failed')

unit.test(all(colnames(W) == rownames(W)),'row and column names of conn mat are same','ERROR with conn mat names')
region.names <- colnames(W) # store connectivity matrix names

####################################################
### Generate ABA-CNDR lists that serve as "keys" ###
####################################################

path.names <- colnames(path.data)[-c(1:2)]
n.regions.CNDR <- length(path.names)
# This key lists the ABA regions associated with each CNDR region
# allowing you to go from a vector of ABA pathology to a vector of CNDR pathology
# The raw key files ('ABAtoCNDR.csv' and 'CNDRtoABA.csv') don't have ipsi and contra for ABA
# since it is assumed that there is perfect L-R symmetry in the ABA connectivity matrix
# so need to loop through and label each ABA region 'i' or 'c' to match the connectivity matrix labels
# I've generated above

# ABA -> CNDR key
# element names: CNDR regions
# element contents: corresponding ABA regions

ABA.to.CNDR.key <- list()
for(CNDR.region in path.names){
  # add i or c to ABA regions
  ABA.region.match <- unique(unlist(strsplit(CNDR.to.ABA$ABA[CNDR.to.ABA$Designation == CNDR.region],','))) # there are duplicates in here so use unique to remove
  hemi <- substr(CNDR.region,1,1)
  ABA.to.CNDR.key[[CNDR.region]] <- paste0(hemi,trimws(ABA.region.match)) # remove any leading or trailing spaces with trimws()
}

# CNDR -> ABA key
# element names: ABA regions
# element contents: corresponding CNDR regions

CNDR.to.ABA.key <- list() 
for(ABA.region in conn.names.ipsi){
  # get the CNDR region(s) that matches the ABA region
  CNDR.region.match <- unlist(strsplit(ABA.to.CNDR$Designation[ABA.to.CNDR$ABA == ABA.region],','))
  if(length(CNDR.region.match)>0){     # if the ABA region has corresponding CNDR regions
    # then store keys for ipsi and contra regions
    CNDR.to.ABA.key[[paste0('i',ABA.region)]] <- paste0('i',trimws(CNDR.region.match))
    CNDR.to.ABA.key[[paste0('c',ABA.region)]] <- paste0('c',trimws(CNDR.region.match))
  } else if(length(CNDR.region.match)==0){
    CNDR.to.ABA.key[[paste0('i',ABA.region)]] <- character(0)
    CNDR.to.ABA.key[[paste0('c',ABA.region)]] <- character(0)
  }
}

# retain indices to reorder like original data variable for plotting on mouse brains
save(path.data, path.names,region.names, CNDR.to.ABA.key, ABA.to.CNDR.key, n.regions.ABA, n.regions.ABA.hemi,n.regions.CNDR,file = paste(savedir,'pathdata.RData',sep=''))
writeMat(paste(savedir,'W.mat',sep=''),W=as.matrix(W))
