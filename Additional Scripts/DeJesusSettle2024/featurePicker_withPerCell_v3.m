% 1/5/22
% Using featurePicker function, generate a structure with protrusion and
% peak stats for all synapses in a given folder of MC synapses
% Fused with features_perCell_v3 on 20230124

% 20230207: Decided w/ MH & AHS to NOT exclude rim relief (c>0) data.

%% Import grids into a cell array 
% This code segment will search a folder containing grids and import
% them into 50x50 data structures and into an array. 

fprintf('Please select folder of interest containing MC grids \n')
mc_file_path = uigetdir('Select Folder of interest');

%Set file paths based on user input
current_dir = pwd;
cd(mc_file_path);
files = dir(mc_file_path);

%Initialize
mc_synapses = {};
names = {};
labels = {};
colors = [];

for file = files'
    if contains(file.name,'_MeanC') % check for only 50x50 grid files
        load(file.name)  %loaded synapse will have variable name 'mc_i'

        mc_synapses = [mc_synapses; mc_i];
        names = [names;file.name]; %store file names for plotting/organization
        
    end
        
end

%% Run through all grids and add stats to nested struct

feature_struct = struct('FileName',{},'Celltype',{},'Protrusion_Stats',{},...
    'Peak_Stats',{},'labeled_Protrusions',{},'labeled_Peaks',{});

for i = 1:length(mc_synapses)
    disp(i)
    % Randomly rotate the synapse before doing the watershedding
    [protrusionStats,peakStats,labeledProtrusions,labeledPeaks] = ...
        featurePicker(imrotate(mc_synapses{i},rand(1)*360,'crop'));
    % Don't randomly rotate
%     [protrusionStats,peakStats,labeledProtrusions,labeledPeaks] = ...
%         featurePicker(mc_synapses{i});    
    feature_struct(i).FileName = names{i};
    feature_struct(i).Protrusion_Stats = protrusionStats;
    feature_struct(i).Peak_Stats = peakStats;
    feature_struct(i).labeled_Protrusions = labeledProtrusions;
    feature_struct(i).labeled_Peaks = labeledPeaks;
    
end

feature_struct = assign_celltype(feature_struct);
feature_struct = sortStruct(feature_struct,'Celltype');

%% Feature-picking for supervised analysis: Protrusions and Peaks
% AHS wrote, MDJ edited 20220718
% V3.0: MDJ: edited 20230124 - Changed peak threshold to >0, and expanded
% calculation of regionprops (peak & protrusion stats)

feature_table = struct2table(feature_struct);

%% Averaging feature data by cell (MDJ)
% Borrowed from MDJ's function for population-averaging in
% analyze_cvurvaturedistancedata_v4.m

features_perCell = struct('FileName',{},'Celltype',{},...
                            'NumProt',[],'NumPeak',[],...
                            'AreaProt_Total',[],'AreaPeak_Total',[],...
                            'AreaProt_Ave',[],'AreaPeak_Ave',[],...
                            'Prot_CurvOverArea',[],'Peak_CurvOverArea',[],...
                            'Prot_CurvMaxAve',[],'Peak_CurvMaxAve',[]);

for idx=1:length(feature_struct)
    
    % Fetch protrusion and peak data
    data_Peak = feature_struct(idx).Peak_Stats;
    data_Prot = feature_struct(idx).Protrusion_Stats;
    
%     % Filter out peripheral rim data because those are different from
%     % peaks: > 80 % of centroid-to-edge ([50 50]), which is ~56.5685 pixel
%     % units
%         centroids_peak = [data_Peak.Centroid];
%         centroids_peak = reshape(centroids_peak,2,size(data_Peak,1))';
%         distances = zeros(1,size(data_Peak,1));
%         for k=1:size(centroids_peak,1)
%             distances(k) = norm([50 50] - centroids_peak(k,:));
%         end
%     % Remove the peak data that are likely to have come from the rim.
%         idx_remove = find(distances>56.568 == 1);
%         data_Peak(idx_remove)=[];
%     % Repeat for prot data
%         centroids_prot = [data_Prot.Centroid];
%         centroids_prot = reshape(centroids_prot,2,size(data_Prot,1))';
%         distances = zeros(1,size(data_Prot,1));
%         for k=1:size(centroids_prot,1)
%             distances(k) = norm([50 50] - centroids_prot(k,:));
%         end
%     % Remove the protrusion data that are likely to have come from the rim.
%         idx_remove = find(distances>56.568 == 1);
%         data_Prot(idx_remove)=[];   
        
    % Fetch and calculate data field by field
    features_perCell(idx).FileName = feature_struct(idx).FileName;
    features_perCell(idx).Celltype = feature_struct(idx).Celltype;
    features_perCell(idx).NumProt = size(data_Prot,1);
    features_perCell(idx).NumPeak = size(data_Peak,1);
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

end
