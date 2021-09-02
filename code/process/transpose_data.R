#################
### Load data ###
#################

rm(list=setdiff(ls(),'params'))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(basedir,'data/')

data <- read.csv(paste0(savedir,'TauPathDataRaw.csv'),stringsAsFactors = F)
data <- t(data)
colnames(data) <- data['Region',]
data <- as.data.frame(data[-which(rownames(data)=='Region'),])
data$Condition <- substr(rownames(data),start=1,stop=3) # extract group label from row names
data$Month <- as.numeric(substr(rownames(data),start=5,stop=5)) # extract MPI
#data$Mouse <- substr(rownames(data),start=8,stop=8) # mouse index

# remove mice that had deeper injections incidentally, per Mike
#mice.remove <- c('G20.3M.5','NTG.6M.1','NTG.6M.4','G20.6M.1','G20.6M.3','NTG.9M.4')
#data <- data[!rownames(data) %in% mice.remove,]

# reorder columns
region.names <- setdiff(colnames(data),c('Condition','Month'))
data <- data[,c('Condition','Month',region.names)]

write.csv(x=data,file = paste0(savedir,'PathData.csv'),row.names = F)
