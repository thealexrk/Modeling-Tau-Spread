# copy over figures for mike -- in order of manuscript

rm(list=ls())
basedir <- '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
setwd(basedir)
injection.site <- c('iDG', 'iCA1', 'iCA3', 'iVISam', 'iRSPagl')
injection.site.label <- paste0(injection.site,collapse='-')
c.max = 0.2
opdir <- paste('TauDiffusion032520_Inject',injection.site.label,'_CMax',c.max,'/',sep='')
fig.dir <- paste0(basedir,'figures_manuscript/')
dir.create(fig.dir,recursive = T)

# Figure 4: bidirectional model fits for NTG, 
grp <- 'NTG'
fig.dir.j <- paste0(fig.dir,'Figure4','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure4') # prefix for figure direction/Figure
desc <- c('script, a: code/diffmodel/optim_bidirectionalspread_onelm_CNDRspace.R and code/diffmodel/optim_bidirectionalspread_onelm_CNDRspace.R',
          'script, b: code/nullmodels/analyzespread_Euclidean_CNDRspace.R and code/nullmodels/plotCNDRspaceEuclideanfit.R',
          'script, d: code/modelcomparison/modelcomparison_traintest.R and code/modelcomparison/plot_modelcomparison_testset.R',
          grp,'(a) A combination of retrograde and anterograde diffusion models explains pathology spread.',
          '(a_labeled) (a), with individual regions labeled.',
          '(b) Predicted vs. actual for optimized diffusion model where the network is the matrix of inverse euclidean distances between the center of mass of each ABA region.',
          '(b_labeled) (b), with individual regions labeled.',
          '(c) alternate seeds fit worse than real seeds for 3-9 MPI',
          '(d) Distributions of model fit in test set using retrograde, anterograde, euclidean, and bidirectional models.',
          '(d, extra) Matrix of fit differences (y-axis minus x-axis) Pairwise one-tailed non-parametric tests computing a p-value for the null hypothesis that \"model on the y-axis fits worse than model on x-axis\"',
          'you can see that all connectivity models beat euclidean distance, then bidirectional > retro > antero.',
          'see csv file for stats')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceFit_bidirectionaladditivemodel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceFit_bidirectionaladditivemodel_label_repel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a_labeled.pdf'))
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'log10predictedpath_bidirectional.csv')
file.copy(from=init,to=paste0(fig.dir.j,'a_log10predicted.csv'))
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean_label_repel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b_labeled.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'SeedSpecificity_OneLM.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'ModelComparisonTestSetPearsonR_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'ModelComparisonTestSetPearsonR_Matrix_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d_diffmatrix.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/',grp,'PValsModelComparisonTestSet.csv')
file.copy(from=init,to=paste0(fig.dir.j,'dStats.csv'))

# Figure 5. Genes and vulnerability
fig.dir.j <- paste0(fig.dir,'Figure5','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure5') # prefix for figure direction/Figure
desc <- c('script, a-: code/diffmodel/plotCNDRspacebidirectionalonelmfit.R',
          'script, c: code/genes_vulnerabiltiy/OteroGarcia_vulnerability_bidirectional.R',
          grp,'using bidirectional one lm vulnerability',
          '(a) hemisphere and (excluding 1 MPI) averaged vulnerability from 4A',
          '(b) hemisphere and time (excluding 1 MPI) averaged vulnerability from 4A vs hemisphere averaged Mapt expression',
          '(c) List of genes whose spatial expression patterns are correlated with NTG vulnerability',
          '(c) Spatial similarity of expression patterns between vulnerability-related genes (FDR corrected p<0.05 cut off for inclusion)',
          '(d) Pearson correlation between all genes and vulnerability, averaged over hemi and time excluding 1 MPI, labeled with pearson r and FDR corrected p-value over ~4200 genes.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'vulnerability_bidirectional_hemiaverage_exclude1 MPI.csv')
file.copy(from=init,to=paste0(fig.dir.j,'a_hemiaveragevulnerabilityexclude1MPI.csv'))
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceHemiAverageVulnerability_Exclude1 MPI_bidirectionalmodel_vsMapt.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/',grp,'Top100GenesVulnerabilitypearson.csv')
file.copy(from=init,to=paste0(fig.dir.j,'cVulnGenes.csv'))
init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/',grp,'SimilarityOfVulnerabilityRelatedGenespearson.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
# init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/',grp,'OteroGarciaGenesFDRSpearman.pdf')
# file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
# init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/',grp,'OteroGarciaGenes.csv')
# file.copy(from=init,to=paste0(fig.dir.j,'c.csv'))
init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/',grp,'VulnerabilityRelatedGenesFDRpearson.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))

# Figure 6 in silico injections
file.copy(from=paste0(opdir,'insilico_injections_onelm'),to=paste0(fig.dir),recursive = T)
file.rename(from=paste0(fig.dir,'insilico_injections_onelm'),to=paste0(fig.dir,'Figure6'))

# figure 8. bootstrap bidirectional models and compare NTG and G20 fits, time constants, anterograde/retrograde betas
fig.dir.j <- paste0(fig.dir,'Figure8','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'Figure8') # prefix for figure direction/Figure
desc <- c('script, a: code/diffmodel/optim_bidirectionalspread_onelm_CNDRspace.R and code/diffmodel/plotCNDRspacebidirectionalonelmfit.R',
          'script, c-e: code/modelcomparison/bootstrap_optimspread_bidirectional.R and code/modelcomparison/G20vsNTG_bootstrap_onelm.R',
          'script, b: code/G20vsNTG/NTGvuln_vs_groupdiff.',
          '(a) In G20 mice, a combination of retrograde and anterograde diffusion models explains pathology spread.',
          '(b) ratio of G20 to NTG regional pathology plotted against hemisphere and time averaged vulnerability, excluding 1 MPI',
          '(c) Distributions of model fit (pearson r) for fitting data to bootstrap samples of mice. NTG and G20 do not differ in model fit (non-parametric, two-tailed test).',
          '(d) Distributions of time constants reveals greater inter-sample variability in retrograde constants compared to anterograde. G20 and NTG do not differ wrt time constants. See Figure*bStats.txt for values. Implies that anterograde may be constant background, but retrograde may hold more importance for explaining individual differences in disease progression and possibly therapeutics as well',
          '(e) Comparing anterograde and retrograde betas between NTG and G20. p-values are bonferroni corrected over all between-group comparisons. anterograde spread importance differs at 6 MPI.',
          '(e, 2) see file Figure8cStats.csv for comparisons between antero and retro within NTG or G20. See critical p-value after bonferroni correction.')#,

write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/G20CNDRSpaceFit_bidirectionaladditivemodel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'G20vsNTG/path_vs_vuln/bidirectional_onelm/',injection.site.label,'/G2019-NTGvsNTGVulnerability_Bidirectional_TimeDependent.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20BootstrapPearson r_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20TimeConstants_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/',grp,'TimeConstantStats.txt')
file.copy(from=init,to=paste0(fig.dir.j,'dStats.txt'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20AnterogradeRetrogradeBetas_Boxplot_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'e.pdf'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGandG20_AnteroVsRetroBetas_Stats.csv')
file.copy(from=init,to=paste0(fig.dir.j,'eStatsAnteroVsRetro.csv'))
init <- paste0(opdir,'modelcomparison/bootstrap/',injection.site.label,'/NTGvsG20_AnteroAndRetroBetas_Stats.csv')
file.copy(from=init,to=paste0(fig.dir.j,'eStatsNTGvsG20.csv'))

# Figure S4: Pathology vs. connectivity

fig.dir.j <- paste0(fig.dir,'FigureS4','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'FigureS4') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c('script, b-c: code/diffmodel/connectivity_vs_path.R',
          grp,'(b-c) Correlation between regional pathology and each region\'s mean anterograde (b) or retrograde (a) connectivity with injection site regions.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'connectivity_vs_path/',injection.site.label,'/',grp,'_AnterogradeConnectivityToInjectionSites_vs_Pathology.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'connectivity_vs_path/',injection.site.label,'/',grp,'_RetrogradeConnectivityToInjectionSites_vs_Pathology.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))

# Figure S5: Pathology vs connectivity for multiple sites

# Figure S6: Diffusion models hemi highlighted, in- out-degree model, conn sim vs alt fit

fig.dir.j <- paste0(fig.dir,'FigureS6','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'FigureS6') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c(grp,
          'script, a: code/diffmodel/plotCNDRspacebidirectionalonelmfit.R',
          'script, b-c: code/nullmodels/optimspread_rewire_onelm_CNDRspace.R and code/nullmodels/plotCNDRspace_rewirefit_onelm.R',
          'script, d: code/nullmodels/seedspec_multi_onelm.R and code/nullmodels/plotseedspec_multi.R',
          '(a) A combination of retrograde and anterograde diffusion models explains pathology spread, with regions colored by their hemisphere.',
          '(b-c) Fits obtained using BIDIRECTIONAL spread along a rewired network that preserves in-degree (b) or out-degree (c), using one lm model.',
          '(d) Alternate seed fit is partially explained by in projection similarity, out projection similarity, distance to injection site,',
          'and distance between alternate seed sites')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/bidirectional_onelm/',injection.site.label,'/',grp,'CNDRSpaceFitHemiColor_bidirectionaladditivemodel.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'nullmodels/rewire_onelm/',injection.site.label,'/',grp,'CNDRSpaceFit_InDegreePreserved_OneLM.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'nullmodels/rewire_onelm/',injection.site.label,'/',grp,'CNDRSpaceFit_OutDegreePreserved_OneLM.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'c.pdf'))
init <- paste0(opdir,'nullmodels/seedspec_multi/',injection.site.label,'/',grp,'RandomSeedFitsVsConnectivity_OneLM.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'d.pdf'))

# Figure S7: consistency of vulnerability across hemispheres and time, genes

fig.dir.j <- paste0(fig.dir,'FigureS7','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'FigureS7') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c('script, a: code/diffmodel/vuln_hemi_time.R',
          'script, b-e: code/genes_vulnerability/genes_vulnerability_bidirectional.R',
          grp,'(a) Spatial similarity of model residuals between each hemisphere-time point combo.',
          '(d) List of genes whose spatial expression patterns are correlated with NTG-G20 pathology',
          '(e) Spatial similarity of expression patterns between NTG-G20 pathology-related genes (FDR corrected p<0.05 cut off for inclusion)',
          '(e) there are ~1200 of these genes so I removed the axis labels but there are two clear clusters')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'diffmodel/vuln_time_hemi/',injection.site.label,'/',grp,'VulnerabilityByTimeAndHemisphere.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/Top100GenesNTGMinusG20PathPearson.csv')
file.copy(from=init,to=paste0(fig.dir.j,'dNTGG20Genes.csv'))
init <- paste0(opdir,'genes_vulnerability/',injection.site.label,'/SimilarityOfNTGMinusG20PathRelatedGenespearson.png')
file.copy(from=init,to=paste0(fig.dir.j,'e.png'))

# Figure S7: consistency of vulnerability across hemispheres and time, genes

fig.dir.j <- paste0(fig.dir,'FigureSX','/') # make figure directory
dir.create(fig.dir.j,recursive = T)
fig.dir.j <- paste0(fig.dir.j,'FigureSX') # prefix for figure direction/Figure
grp <- 'NTG'
desc <- c('script, a: code/nullmodels/analyzespread_Euclidean_CNDRspace.R and code/nullmodels/plotCNDRspaceEuclideanfit.R',
          'script, b-e: code/modelcomparison/modelcomparison_traintest_exclinj.R and code/modelcomparison/plot_modelcomparison_testset_exclinj.R',
          grp,'(a) Predicted vs. actual for optimized diffusion model where the network is the matrix of inverse euclidean distances between the center of mass of each ABA region, excluding injection sites.',
          '(b) Distributions of model fit in test set using euclidean model, euclidean model excluding injection sites, and bidirectional model.')
write.table(x=desc,file = paste0(fig.dir.j,'.txt'),sep = '\n',row.names = F,col.names = F)
init <- paste0(opdir,'nullmodels/euclidean/',injection.site.label,'/',grp,'CNDRSpaceFit_Euclidean_ExcludeInjectionSites.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'a.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/euclidean_exclinj/',grp,'EuclideanExclInjVsAllPoints_TestSetPearsonR_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/euclidean_exclinj/',grp,'EuclideanExclInjVsAllPoints_PValMatrix_TestSetPearsonR_CNDRSpace.pdf')
file.copy(from=init,to=paste0(fig.dir.j,'b_diffmatrix.pdf'))
init <- paste0(opdir,'modelcomparison/traintest/',injection.site.label,'/euclidean_exclinj/',grp,'PValsEuclideanExclInjTestSet.csv')
file.copy(from=init,to=paste0(fig.dir.j,'bStats.csv'))

#########



