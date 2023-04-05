function cellStats = getFeatureStatsperCell(protrusionStats,peakStats,mc_i,scale)
% Given featureStats output from featurePicker, generate a per-cell data
% structure that summarizes the propertise of the features
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023
%
% Arguments:
%       protrusionStats: a MATLAB struct, output from featurePicker.m,
%       defining the negative curvature zones of interest
%       peakStats: a MATLAB struct, output from featurePicker.m,
%       defining the positive curvature zones of interest
%       
% Output:
%
features_perCell = struct('NumProt',[],'NumPeak',[],...
                            'AreaProt_Total',[],'AreaPeak_Total',[],...
                            'AreaProt_Ave',[],'AreaPeak_Ave',[],...
                            'Prot_CurvOverArea',[],'Peak_CurvOverArea',[],...
                            'Prot_CurvMaxAve',[],'Peak_CurvMaxAve',[]);
idx=1;

% Fetch protrusion and peak data
data_Peak = scaleFeatures(peakStats,mc_i,scale);
data_Prot = scaleFeatures(protrusionStats,mc_i,scale);

% Removed this section from final version.
% % Filter out peripheral rim data because those are different from
% % peaks: > 80 % of centroid-to-edge
% grid_center = [scale/2 scale/2];
% rim_threshold = 0.8*scale;
% 
% centroids_peak = [data_Peak.Centroid];
% centroids_peak = reshape(centroids_peak,2,size(data_Peak,1))';
% distances = zeros(1,size(data_Peak,1));
% for k=1:size(centroids_peak,1)
%     distances(k) = norm(grid_center - centroids_peak(k,:));
% end
% % Remove the peak data that are likely to have come from the rim.
% idx_remove = find(distances>rim_threshold == 1);
% data_Peak(idx_remove)=[];
% % Repeat for prot data
% centroids_prot = [data_Prot.Centroid];
% centroids_prot = reshape(centroids_prot,2,size(data_Prot,1))';
% distances = zeros(1,size(data_Prot,1));
% for k=1:size(centroids_prot,1)
%     distances(k) = norm(grid_center - centroids_prot(k,:));
% end
% % Remove the protrusion data that are likely to have come from the rim.
% idx_remove = find(distances>rim_threshold == 1);
% data_Prot(idx_remove)=[];   

% Fetch and calculate data field by field
if ~isempty(data_Prot)
    features_perCell(idx).NumProt = size(data_Prot,1);
else
    features_perCell(idx).NumProt = 0;
end

if ~isempty(data_Peak)
    features_perCell(idx).NumPeak = size(data_Peak,1);
else
    features_perCell(idx).NumPeak = 0;
end

areas_prot = [data_Prot.Area];
    features_perCell(idx).AreaProt_Total = sum(areas_prot);
    features_perCell(idx).AreaProt_Ave = mean(areas_prot);
areas_peak = [data_Peak.Area];
    features_perCell(idx).AreaPeak_Total = sum(areas_peak);
    features_perCell(idx).AreaPeak_Ave = mean(areas_peak);
MFIcurvatures_prot = [data_Prot.MeanIntensity];
    features_perCell(idx).Prot_CurvOverArea = sum(areas_prot.*MFIcurvatures_prot);
MFIcurvatures_peak = [data_Peak.MeanIntensity];
    features_perCell(idx).Peak_CurvOverArea = sum(areas_peak.*MFIcurvatures_peak);
Maxcurvatures_prot = [data_Prot.MinIntensity];    
    features_perCell(idx).Prot_CurvMaxAve = mean(Maxcurvatures_prot);
Maxcurvatures_peak = [data_Peak.MaxIntensity];
    features_perCell(idx).Peak_CurvMaxAve = mean(Maxcurvatures_peak);

cellStats = features_perCell;

end

