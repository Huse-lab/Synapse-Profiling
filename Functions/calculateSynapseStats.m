function synapseResults = calculateSynapseStats(MPStats,roidata)
% Unifying function to generate a MATLAB struct containing all the
% statistics and information analyzed as part of de Jesus and Settle et al.
% This runs on one MPStats element (can be one part of the array) and one
% ROI data structure 
%
% Arguments
%      - MPStats: a 1x1 MPStats/MPRender struct, output from Preprocessing
%      - roidata: a 1x1 ROIdata structure, output from selectROI_batch.m
% Outputs
%      - see READMe for complete explanation of the synapseResults
%      structure


%check that arguments are good
if size(MPStats,1) > 1 
    error(['Error: calculateSynapseStats can only handle a 1x1 struct....' ...
        'To calculate in batch, loop through the MPStats structure one at a time']) 
elseif size(roidata,1) > 1
    error(['Error: calculateSynapseStats can only handle a 1x1 struct....' ...
        'To calculate in batch, loop through the ROI structure one at a time']) 
end


%Run the barrage of analysis functions
[mc_i,gc_i,stain_i,scale_um] = convertROI_toGrid(MPStats,'ROI',roidata,50);
LPA_output = calculateBulkDeformation(MPStats,roidata,20);
Z_coefficients = calculateZernikeCoefficients(mc_i,15);
RDZ_coefficients = deRotateZernikeCoefficients(Z_coefficients);
[Zernike_Complexity, fit_trace] = calculateZernikeComplexity(mc_i,0.8);
[~,edgecoor_sph,mc,gc] = calculateCurvatures(MPStats);
curvatureStats = calculateCurvatureDistribution(MPStats,roidata);
mcMatrix = calculate_DistanceCurvatureMatrix(curvatureStats,10,(-1:0.1:1),'mc');
gcMatrix = calculate_DistanceCurvatureMatrix(curvatureStats,10,(-0.3:0.03:0.3),'gc');

[protrusionStats,peakStats,labeledProtrusions,labeledPeaks] = featurePicker_noFilter(mc_i);
featureStats = getFeatureStatsperCell(protrusionStats,peakStats,mc_i,scale_um);

synapseResults = struct();
synapseResults.FileName = MPStats.FileName;
synapseResults.MeanCurvatureGrid = mc_i;
synapseResults.GaussianCurvatureGrid = gc_i;
synapseResults.StainIntensityGrid = stain_i;
synapseResults.gridScale_um = scale_um; %width of grid in µm
synapseResults.ZernikeModalCoefficients = Z_coefficients;
synapseResults.ZernikeRDCoefficients = RDZ_coefficients;
synapseResults.ZernikePatternComplexity = Zernike_Complexity;
synapseResults.ZernikeComplexityTrace = fit_trace;
synapseResults.RadialCurvatureData = curvatureStats;
synapseResults.MeanCurvatureRadialProfile = mcMatrix;
synapseResults.GaussianCurvatureRadialProfile = gcMatrix;
synapseResults.BulkDeformationResults = LPA_output;
synapseResults.ProtrusonStats = protrusionStats;
synapseResults.PeakStats = peakStats;
synapseResults.FeatureSummary = featureStats;


end