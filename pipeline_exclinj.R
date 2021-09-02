rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
setwd(basedir)
params <- list(basedir=basedir,      
               matlab.path='/Applications/MATLAB_R2019a.app/bin/matlab',
               grps = c('NTG','G20'),
               #injection.site=c('iDG','iVISam'), # define injection sites (using ABA nomenclature)
               injection.site = c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl'), # injection sites (using ABA nomenclature)
               injection.site.CNDR = c(iDG='iDG', iCA1='iCA1', iCA3='iCA3', iVISam='iVISa', iRSPagl='iRSP'), # convert ABA injection sites to CNDR
               tps=c(1,3,6,9), # define time points post injection
               c.min = 1e-5, # empirically conservative minimum for time constant
               c.max = 0.2, # empirically conservative maximum for time constant
               c.n = 100) # number of time constants to try
params$excl.inj.CNDR <- unname(params$injection.site.CNDR[params$injection.site]) # exclude injection sites from fit calculation; if NULL all injection sites are included
source('code/misc/miscfxns.R')
params$source.save <- source.save
params$opdir <- paste('TauDiffusion100120_Inject',paste0(params$injection.site,collapse='-'),'_CMax',params$c.max,'/',sep='')
dir.create(params$opdir,recursive = T)

#################################################
### Load packages & create output directories ###
#################################################

source('code/misc/packages.R')
# start up for many scripts
injection.site <- params$injection.site; grp <- 'NTG'
##############################
### Process pathology data ###
##############################

source('code/process/transpose_data.R')
source('code/process/process.R')
#source('code/process/testCNDRtoABAmatching.R')

########################################
### Process ABA gene expression data ###
########################################

source('code/aba/process_ontology.R')
source('code/aba/download_gene_expression.R') # specify variable gois in this file to get expression for whatever genes you want

###################################################
### Get euclidean distances between ABA regions ###
###################################################

source('code/aba/atlas_structures.R') # required for most scripts in code/nullmodels/

#######################
### Diffusion model ###
#######################

# retrograde model with additive Mapt expression
injection.site <- params$injection.site
grp <- 'NTG'
source('code/diffmodel/analyze_retrogradespread_CNDRspace.R')
goi <- 'Mapt'
probe <- 'RP_071204_01_D02'
source('code/diffmodel/plotCNDRspacefit.R')

# bidirectional diffusion model with different a single linear model weighting anterograde and retrograde
# Manuscript: Figure 4A (NTG) model fit by month, 5B (NTG) residuals vs. MAPT exp
for(grp in params$grps){
  injection.sites <- c(list(params$injection.site))
  for(injection.site in injection.sites){
    source('code/diffmodel/optim_bidirectionalspread_onelm_CNDRspace.R')
    goi <- 'Mapt'
    source('code/diffmodel/plotCNDRspacebidirectionalonelmfit.R')
  }
}
################################
### Quality control analyses ###
################################

# these are very time consuming

for(grp in params$grps){
  # seed specificity
  #source('code/diffmodel/seedspec.R')
  #source('code/diffmodel/plotseedspec.R')
  
  # validate time constant out of sample
  #source('code/diffmodel/traintest_bidirectional_CNDRspace.R')
  source('code/diffmodel/plot_traintest_bidirectional.R')
}

# bootstrap time constants and fits
injection.sites <- list(params$injection.site)
for(injection.site in injection.sites){
  for(grp in params$grps){
    source('code/modelcomparison/bootstrap_optimspread_bidirectional.R')
  }
  grp <-'NTG'
  #source('code/modelcomparison/modelcomparison_traintest.R')
}
source('code/modelcomparison/plot_modelcomparison_testset.R')
source('code/modelcomparison/plot_modelcomparison_testset_exclinj.R')
###############################
### Genes and vulnerability ###
###############################

injection.site <- params$injection.site; grp <- 'NTG'
source('code/genes_vulnerability/genes_vulnerability_bidirectional.R')

###########################
### Network null models ###
###########################

# use each seed site separately, all together, and then do entire hippocampus
#injection.sites <- c(as.list(params$injection.site),list(params$injection.site),list(c('iDG','iCA1','iCA3')))
injection.sites <- c(list(params$injection.site))
# retrograde model with additive Mapt expression
for(injection.site in injection.sites){
  for(grp in params$grps){
    source('code/nullmodels/analyzespread_Euclidean_CNDRspace.R')
    goi <- 'Mapt'
    probe <- 'RP_071204_01_D02'
    source('code/nullmodels/plotCNDRspaceEuclideanfit.R')
  }
}

# run matlab scripts to generate distance matrix and rewired connectivity matrix
mat.savedir <- paste(params$basedir,params$opdir,'processed',sep='')
mat.cmd <- paste(params$matlab.path,' -nojvm -r \"cd(\'',params$basedir,'code/nullmodels/\'); homedir = \'',params$basedir,'\'; ',
                 'savedir = \'', mat.savedir,'\'; run(\'makenullconn.m\'); exit;\"',sep='')
system(mat.cmd)

grp <- 'NTG'
injection.site <- params$injection.site
source('code/nullmodels/optimspread_rewire_onelm_CNDRspace.R')
source('code/nullmodels/plotCNDRspace_rewirefit_onelm.R')
############################
### G20 vs. NTG analyses ###
############################

# use each seed site separately, all together, and then do entire hippocampus
injection.sites <- c(as.list(params$injection.site),list(params$injection.site),list(c('iDG','iCA1','iCA3')))
for(injection.site in injection.sites){
  source('code/G20vsNTG/NTGvuln_vs_groupdiff.R')
}

#####################
### miscellaneous ###
#####################

# use each seed site separately, all together, and then do entire hippocampus
injection.sites <- c(as.list(params$injection.site),list(params$injection.site),list(c('iDG','iCA1','iCA3')))
for(injection.site in injection.sites){
  for(grp in params$grps){
    #source('code/diffmodel/connectivity_vs_path.R')
    source('code/diffmodel/vuln_hemi_time.R')
  }
}

# in silico injections with bidirectional model time constants and regression weights
injection.sites <- c('iENTl', 'iSNc', 'iCP','iXII','iNTS'); grp <- 'NTG'
for(injection.site in injection.sites){source('code/diffmodel/insilico_inject.R')}

# in silico injections with bidirectional model that only has one set of weights and time constants
injection.sites <- c('iENTl', 'iSNc', 'iCP','iXII','iNTS'); grp <- 'NTG'
for(injection.site in injection.sites){source('code/diffmodel/insilico_inject_onelm.R')}
