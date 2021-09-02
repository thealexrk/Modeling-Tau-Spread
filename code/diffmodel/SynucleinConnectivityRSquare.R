#################
### Load data ###
#################

rm(list=setdiff(ls(),c('params','grp','goi','probe')))
print(paste('Running:',grp,goi,probe))
basedir <- params$basedir
setwd(basedir)
savedir <- paste0(params$opdir,'diffmodel/')
dir.create(savedir,recursive=T)

source('code/misc/fitfxns.R')
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names
gene.exp <- read.csv('data/aba/expression/GeneExpressionABA.csv') # load gene expression
gene.exp <- as.matrix(gene.exp[,paste0(goi,'_',probe)])
rownames(gene.exp) <- region.names

# get mean of each month
tp <- 1
Mice <- path.data[path.data$Condition == grp,-1]
Grp.mean <- colMeans(Mice,na.rm = T)
log.path <- log(Grp.mean,base=10)
load(file = paste0(savedir,grp,'CNDRSpaceFit_data.RData'))
Xt.base <- df$pred
load(file = paste0(savedir,goi,probe,'_Weighted/',grp,goi,probe,'Weighted_CNDRSpaceFit_data.RData'))
Xt.gene <- df$pred

######################################
### Connectivity-only "base" model ###
######################################

# make data frame of log path, log path prediction, and synuclein
df <- data.frame(path = log.path, Xt = Xt.base, Gene = map.ABA.to.CNDR(gene.exp,path.names,ABA.to.CNDR.key))
# exclude regions with 0 pathology at each time point for purposes of computing fit
mask <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
print(paste(sum(mask),'regions')) # number of regions left after exclusion
df <- df[mask,] 
m.conn <- lm(path~Xt,data=df)
print(paste('diffusion R^2 =',summary(m.conn)$r.squared))
m.gene <- lm(path~Gene,data=df)
print(paste('synuclein R^2 =',summary(m.gene)$r.squared))

full.mdl <- lm(path~Xt + Gene,data=df)
print(paste('max possible R^2 =',summary(full.mdl)$r.squared))
print(paste('connectivity-independent variance in pathology explained by synuclein:',cor(residuals(m.conn),df$Gene)^2))
print(paste('synuclein-independent variance in pathology explained by connectivity:',cor(residuals(m.gene),df$Xt)^2))

################################
### Synuclein-weighted model ###
################################

# make data frame of log path, log path prediction, and synuclein
df <- data.frame(path = log.path, Xt = Xt.gene, Gene = map.ABA.to.CNDR(gene.exp,path.names,ABA.to.CNDR.key))
# exclude regions with 0 pathology at each time point for purposes of computing fit
mask <- df$path != -Inf & df$Xt != -Inf & !is.na(df$Xt)
print(paste(sum(mask),'regions')) # number of regions left after exclusion
df <- df[mask,] 
m.conn <- lm(path~Xt,data=df)
print(paste('syn-weighted diffusion R^2 =',summary(m.conn)$r.squared))
m.gene <- lm(path~Gene,data=df)
print(paste('synuclein R^2 =',summary(m.gene)$r.squared))

full.mdl <- lm(path~Xt + Gene,data=df)
print(paste('max possible R^2 =',summary(full.mdl)$r.squared))
print(paste('syn-weighted connectivity-independent variance in pathology explained by synuclein:',cor(residuals(m.conn),df$Gene)^2))
print(paste('synuclein-independent variance in pathology explained by syn-weighted connectivity:',cor(residuals(m.gene),df$Xt)^2))
