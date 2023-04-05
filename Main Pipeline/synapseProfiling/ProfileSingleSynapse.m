% Script to load a single pair of MPRender and ROI data files and
% calculate all relevant stats. See manuscript for general description of
% methodology and the "calculateSynapseResults.m" function for specific
% outputs. 
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023

%% Load MPStats and ROI, then run 
[MPfile, MPfilepath] = uigetfile('.mat',"Choose the MPStats/MPRender File");
[ROIfile, ROIfilepath] = uigetfile('.mat',"Choose the MPStats/MPRender File");

MPhandle = load(fullfile(MPfilepath,MPfile));
ROIhandle = load(fullfile(ROIfilepath,ROIfile));

MPStats = MPhandle.MPStats;
%to smooth any inconsistencies in variable naming
roidata_Field = fieldnames(ROIhandle);
roidata = ROIhandle.(roidata_Field{1});


synapseResults = calculateSynapseStats(MPStats,roidata);
generateSynapseResultsPage(MPStats,roidata,synapseResults);
outpng = strrep(MPStats.FileName,'.tif',' - ResultsOverview.png');
saveas(gcf,outpng)
close
outmat = strrep(MPStats.FileName,'.tif',' - ResultStuct.mat');
save(outmat,'synapseResults');