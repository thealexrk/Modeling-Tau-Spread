% this code is borrowed from http://help.brain-map.org/display/mousebrain/API#API-DownloadAtlas
% and modified by EJC on 4-11-20
clear all; close all; clc
basedir = '~/Dropbox/Neurodegeneration/TauSpread/tau-spread/';
cd(basedir)
%%
% Download and unzip the atlasVolume and annotation zip files
 
% 25 micron volume size
size = [528 320 456];
% VOL = 3-D matrix of atlas Nissl volume
fid = fopen('data/aba/atlas/atlasVolume/atlasVolume.raw', 'r', 'l' );
VOL = fread( fid, prod(size), 'uint8' );
fclose( fid );
VOL = reshape(VOL,size);
% ANO = 3-D matrix of annotation labels
fid = fopen('data/aba/atlas/P56_Mouse_annotation/annotation.raw', 'r', 'l' );
ANO = fread( fid, prod(size), 'uint32' );
fclose( fid );
ANO = reshape(ANO,size);
ANOunique = unique(ANO);

% 200 micron volume size
sizeGrid = [67 41 58];
% ANOGD = 3-D matrix of grid-level annotation labels
fid = fopen('data/aba/atlas/P56_Mouse_gridAnnotation/gridAnnotation.raw', 'r', 'l' );
ANOGD = fread( fid, prod(sizeGrid), 'uint32' );
fclose( fid );
ANOGD = reshape(ANOGD,sizeGrid);
ANOGDunique = unique(ANOGD);

% save data to process in R
save('data/aba/atlas/atlasABA.mat','VOL','ANO','ANOunique','ANOGD','ANOGDunique','size','sizeGrid');
%% plot
% Display one coronal section
figure;imagesc(squeeze(VOL(264,:,:)));colormap(gray);
figure;imagesc(squeeze(ANO(264,:,:)));colormap(lines); colorbar
 
% Display one sagittal section
figure;imagesc(squeeze(ANO(:,:,220)));colormap(lines);
figure;imagesc(squeeze(VOL(:,:,220)));colormap(gray);

% Display one coronal and one sagittal section
figure;imagesc(squeeze(ANOGD(34,:,:)));colormap(lines);
figure;imagesc(squeeze(ANOGD(:,:,28)));colormap(lines);