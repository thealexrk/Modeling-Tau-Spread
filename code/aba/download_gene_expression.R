rm(list=setdiff(ls(),c('params')))
basedir <- params$basedir # make this whatever you like, end with /
setwd(basedir)
source('code/aba/aba_fxns.R')
source('code/misc/fitfxns.R')
abadir <- 'data/aba/'
dir.create(paste0(abadir,'expression'),recursive = T)
load(file=paste0(abadir,'ontologies/keys.RData'))
load(paste(params$opdir,'processed/pathdata.RData',sep=''))  # load path data and ROI names

dfile <- paste0(basedir,'data/aba/mouse_expression_data_sets.csv')
if(!file.exists(dfile)){
  download.file(url = 'http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/mouse_expression_data_sets.csv',
  destfile = dfile)}
datasets <- read.csv('data/aba/mouse_expression_data_sets.csv')
datasets <- datasets[datasets$plane_of_section == 'coronal',] # only look at coronal sections
# see http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/Accessing_October_2014_expression_data.pdf
# obtained from http://download.alleninstitute.org/informatics-archive/october-2014/mouse_expression/mouse_expression_data_sets.csv

## general framework for downloading data ##
gois <- list(Tau=list(goi='Mapt',probe='RP_071204_01_D02'))
df.expression <- ontology.get.expression.df(region.names,region.names.id.key,gois,datasets,download=TRUE)
write.csv(x = df.expression,file = paste0(abadir,'expression/GeneExpressionABA.csv'),row.names = FALSE)

## gene expression used for this project ##
Tau <- map.ABA.to.CNDR(df.expression,path.names,ABA.to.CNDR.key)
write.csv(x = Tau,file = paste0(abadir,'expression/',gois$Tau$goi,gois$Tau$probe,'ExpressionCNDR.csv'))
