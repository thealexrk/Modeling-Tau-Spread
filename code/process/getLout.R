rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste(params$opdir,'processed/',sep='')
dir.create(savedir,recursive=T)

load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
load(paste(params$opdir,'processed/Snca.RData',sep='')) # load Snca expression

# 3-1-19 Some Synuclein expression values are negative... shifting to 1 as the minimum
# so that if there's no synuclein in a region you don't change that row
# when you weight the matrix

Synuclein <- Synuclein - min(Synuclein,na.rm=T) + 1 # set new minimum to 1

W <- readMat(paste(params$opdir,'processed/W.mat',sep=''))$W
W <- W * !diag(n.regions) # get rid of diagonal
W <- diag(as.numeric(Synuclein)) %*% W
#W <- W / (max(Re(eigen(W)$values))) # scale to max eigenvalue

# Where i is row element and j is column element
# Wij is a connection from region i to region j
# convention is the opposite, so without transposing W
# I am capturing "retrograde" connectivity

in.deg <- colSums(W)
out.deg <- rowSums(W)
L.out <- diag(x = out.deg) - W # outdegree laplacian
save(L.out,file = paste(params$opdir,'processed/Lout.RData',sep=''))
