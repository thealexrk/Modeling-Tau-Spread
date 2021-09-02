# load atlas data
source('code/misc/plottingfxns.R') # I'm going to use my imagesc function here to visualize images
basedir <- '/Users/Eli/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
abadir <- paste0(basedir,'data/aba/')
load(file=paste0(abadir,'ontologies/keys.RData'))
atlas.fname <- paste0(abadir,'atlas/mask_aba.mat') # output filename for mask
pymask <- readMat(atlas.fname)

# extract variables from .mat file
atlas <- pymask$mask # the atlas is a 3D matrix whose elements are either 0 for non-brain voxels, or an integer that identifies a region
inds <- pymask$region.inds # all of the integer indices that identify regions. here 1:213 because they treat every region as hemispherically symmetric

# load region names
load(file=paste0(abadir,'ontologies/keys.RData'))
region.names <- names(region.names.id.key)


# example of general matrix manipulation

X <- matrix(rnorm(100),nrow=20,ncol=5) # matrix of numbers normally distributed around 0
X # check it out
mask <- X < 0# threshold at 0 by creating a logical mask (X < 0)
mask # check it out
X[mask] <- 0 
X # check it out

# alternatively instead of using masks we can use integer indices to replace all the 0's with 3's
inds <- which(X==0,arr.ind = T)
for(j in 1:nrow(inds)){
  X[inds[j,'row'],inds[j,'col']] <- 3  
}
X  # check it out

# okay now on to the atlas

# check dimensions of atlas
print(dim(atlas))
#image.upright <- function(x) image(t(x[nrow(x):1,])) # can also use R's native image function but it flips images by default
#image.upright(atlas[260,,])
# visualize a coronal slice in middle of brain
imagesc(atlas[260,,],noticks = T,cmap = 'Blues') 
# ^-- here the color axis corresponds to these region integer indexes, rather than any data.
# if we wanted to highlight region number 50 and fill it with data, we can do this
idx.plot <- 50 # let's plot region 50
region.names[idx.plot] # turns out that's EPv
atlas.slice <- atlas[260,,] # just take that coronal slice and store in a new matrix
atlas.plot <- atlas.slice==idx.plot # plot a mask that just highlights that region
imagesc(atlas.plot,noticks = T,cmap='Blues') # now all we see is EPv
# hold on this for now

# currently atlas doesn't split left and right hemisphere, so shift region integer identifiers on one hemisphere
atlas.hemi.mask <- atlas*0 # make a new matrix of all 0's 
atlas.hemi.mask[,,(dim(atlas)[3]/2):(dim(atlas)[3])] <- TRUE # set one hemisphere to TRUE and the rest to FALSE 
for(idx in inds){
  print(idx)
  # for each region (which contains both hemispheres), split into a left and right hemi region
  # by adding number of regions in hemisphere to each index
  # atlas == idx makes a logical mask of dimensions same as atlas with values TRUE where the atlas elements equal one regions index
  # atlas.hemi.mask is a static logical mask that highlights one hemisphere at a time
  atlas[atlas == idx & atlas.hemi.mask] <- idx + length(inds) 
}
#^-- this takes a while

region.data <- runif(length(region.names)) # simulate region data uniform from 0 to 1
# populate regions in atlas slice above with data
for(idx in 1:length(region.names)){
  atlas.slice[atlas.slice==idx] <- region.data[idx]
}
imagesc(atlas.slice,noticks = T,cmap='Blues')
