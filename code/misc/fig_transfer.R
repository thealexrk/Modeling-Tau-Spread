# copy over figures for mike

rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
setwd(basedir)
injection.site <- c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl')
injection.site.label <- paste0(injection.site,collapse='-')
c.max = 0.2
opdir <- paste('TauDiffusion032520_Inject',injection.site.label,'_CMax',c.max,'/',sep='')
fig.dir <- paste0(basedir,'figures/')
dir.create(fig.dir,recursive = T)

# Figure 1&2: bidirectional model fits for NTG and G20
f.defs <- list(list(grp='NTG',f='Figure1'),list(grp='G20',f='Figure2'))
for(f.def in f.defs){
  f <- f.def$f
  grp <- f.def$grp
  fig.dir.j <- paste0(fig.dir,f,'/') # make figure directory
  dir.create(fig.dir.j,recursive = T)
  fig.dir.j <- paste0(fig.dir.j,f) # prefix for figure direction/Figure
  desc <- c('script: code/diffmodel/plotCNDRspacebidirectionalonelmfit.R',
            grp,'(a) A combination of retrograde and anterograde diffusion models explains pathology spread.',
            '(b) The same as (a), coloring each point by hemisphere relative to injection site.',
            '(c) hemisphere and time (excluding 1 MPI) averaged vulnerability from (a) vs hemisphere averaged Mapt expression')
  write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
  init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceFit_bidirectionaladditivemodel.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceFitHemiColor_bidirectionaladditivemodel.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceHemiAverageVulnerability_Exclude1 MPI_bidirectionalmodel_vsMapt.pdf')
  file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
  init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'vulnerability_bidirectional_hemiaverage_exclude1 MPI.csv')
  file.copy(from=init,to=paste0(fig.dir.j,'a_hemiaveragevulnerabilityexclude1MPI.csv'))
  init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'log10predictedpath_bidirectional.csv')
  file.copy(from=init,to=paste0(fig.dir.j,'a_log10predicted.csv'))
}

# Figure 3. Retrograde alone: CNDRSpaceCSweepByTimePoint_basemodel.pdf
# shows that average time constant is a good approximation of "real" time constant, 1 MPI not biasing significantly
fig.dir.j <- paste0(fig.dir,'Figure3','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure3') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) For retrograde model, average time constant is a good approximation of "real" time constant, 1 MPI not biasing significantly. vertical line is group average time constant.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/retrograde/',injection.site.label,'/',grp,'CNDRSpaceCSweepByTimePoint_basemodel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))

# Figure 4. Null model of euclidean distance as network
fig.dir.j <- paste0(fig.dir,'Figure4','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure4') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) Predicted vs. actual for optimized diffusion model where the network is the matrix of inverse euclidean distances between the center of mass of each ABA region.',
          '(b) Same plot as (a), only I have excluded the outliers where predicted path > 95th %ile. Empirically this always includes the injection sites')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean_ExcludeInjectionSites.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))

# in silico injections
file.copy(from=paste0(opdir,'insilico_injections_onelm'),to=paste0(fig.dir),recursive = T)

# Figure 5. Relationship between vulnerability and pathology in NTG-G20, and consistency of vulnerability across hemispheres and time.
fig.dir.j <- paste0(fig.dir,'Figure5','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure5') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c('script, a-b: code/G20vsNTG/NTGvuln_vs_groupdiff.',
          'script, c-d: code/diffmodel/vuln_hemi_time.R',
          'All using bidirectional one lm model',
          '(a) ratio of G20 to NTG regional pathology plotted against time-point specific vulnerability',
          '(b) ratio of G20 to NTG regional pathology plotted against hemisphere and time averaged vulnerability, excluding 1 MPI',
          grp,'(c) Vulnerability values compared between hemispheres at each time point using wilcox rank sum. p-vals not adjusted for MC.',
          '(d) Spatial similarity of model residuals between each hemisphere-time point combo.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'G20vsNTG/path_vs_vuln/bidirectional_onelm/',injection.site.label,'/G2019-NTGvsNTGVulnerability_Bidirectional_HemiTimeAverage_exclude1 MPI.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'G20vsNTG/path_vs_vuln/bidirectional_onelm/',injection.site.label,'/G2019-NTGvsNTGVulnerability_Bidirectional_TimeDependent.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'diffmodel/vuln_time_hemi/',injection.site.label,'/',grp,'CompareVulnerabilityBetweenHemispheres.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'diffmodel/vuln_time_hemi/',injection.site.label,'/',grp,'VulnerabilityByTimeAndHemisphere.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))

# alternate seed sites
fig.dir.j <- paste0(fig.dir,'Figure6','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure6') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a) distance between random seed sites is within 10% of distance between tested injection sites',
          '(b) alternate seeds fit worse than real seeds for 3-9 MPI',
          '(c) Alternate seed fit is partially explained by in projection similarity, out projection similarity, distance to injection site,',
          'and distance between alternate seed sites')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/AlternateSeedSiteDistances.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'SeedSpecificity.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'RandomSeedFitsVsConnectivity.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))

# Figure 7. comparing retrograde, anterograde, bidirectional models in unseen test-set data
fig.dir.j <- paste0(fig.dir,'Figure7','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure7') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'scripts: code/modelcomparison/modelcomparison_traintest.R and code/modelcomparison/plot_modelcomparison_testset.R',
          '(a) Distributions of model fit in test set using retrograde, anterograde, euclidean, and bidirectional models.',
          '(b) Matrix of fit differences (y-axis minus x-axis) Pairwise one-tailed non-parametric tests computing a p-value for the null hypothesis that \"model on the y-axis fits worse than model on x-axis\"',
          'you can see that all connectivity models beat euclidean distance, then bidirectional > retro > antero.')

write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'ModelComparisonTestSetPearsonR_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'ModelComparisonTestSetPearsonR_Matrix_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'PValsModelComparisonTestSet.csv')
file.copy(from=init,to=paste0(fig.dir.j,'pvals.csv'))

# figure 8. bootstrap bidirectional models and compare NTG and G20 fits, time constants, anterograde/retrograde betas
fig.dir.j <- paste0(fig.dir,'Figure8','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure8') # prefix for figure direction/Figure
desc <- c('code: code/modelcomparison/bootstrap_optimspread_bidirectional.R and code/modelcomparison/G20vsNTG_bootstrap_onelm.R',
          '(a) Distributions of model fit (pearson r) for fitting data to bootstrap samples of mice. NTG and G20 do not differ in model fit (non-parametric, two-tailed test).',
          '(b) Distributions of time constants reveals greater inter-sample variability in retrograde constants compared to anterograde. G20 and NTG do not differ wrt time constants. See Figure*bStats.txt for values. Implies that anterograde may be constant background, but retrograde may hold more importance for explaining individual differences in disease progression and possibly therapeutics as well',
          '(c) Comparing anterograde and retrograde betas between NTG and G20. p-values are bonferroni corrected over all between-group comparisons. anterograde spread importance differs at 6 MPI.',
          '(c, 2) see file Figure8cStats.csv for comparisons between antero and retro within NTG or G20. See critical p-value after bonferroni correction.')#,
          #'(d) non-parametric comparison of vulnerability values between NTG and G20 reveals significant differences in vulnerability between the groups.',
          #'(d, 2) see Figure8d_pvals.txt for definition of asterisks. p-vals are uncorrected b/c we need a ton of bootstraps to have numerical precision for the new p-crit after bonferroni correction.')

write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20BootstrapPearson r_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20TimeConstants_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/',grp,'TimeConstantStats.txt')
file.copy(from=init,to=paste0(fig.dir.j,'bStats.txt'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20AnterogradeRetrogradeBetas_Boxplot_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGandG20_AnteroVsRetroBetas_Stats.csv')
file.copy(from=init,to=paste0(fig.dir.j,'cStatsAnteroVsRetro.csv'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20_AnteroAndRetroBetas_Stats.csv')
file.copy(from=init,to=paste0(fig.dir.j,'cStatsNTGvsG20.csv'))
# init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20BootstrapVulnerability_CNDRSpace.pdf')
# file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))
# init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20BootstrapVulnerability_CNDRSpace_pkey.txt')
# file.copy(from=init,to=paste0(fig.dir.j,'d_pvals.txt'))

# Figure 9. degree-preserving null models
fig.dir.j <- paste0(fig.dir,'Figure9','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure9') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,'(a-b) Fits obtained using retrograde spread along a rewired network that preserves in-degree (a) or out-degree (b).')

write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/rewire/',injection.site.label,'/',grp,'CNDRSpaceFit_InDegreePreserved.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/rewire/',injection.site.label,'/',grp,'CNDRSpaceFit_OutDegreePreserved.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))


