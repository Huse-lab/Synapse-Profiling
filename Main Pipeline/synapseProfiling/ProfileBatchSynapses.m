% Batch script to load a set of MPRender and ROI data files and
% calculate all relevant stats. See manuscript for general description of
% methodology and the "calculateSynapseResults.m" function for specific
% outputs. 
%
% All MPRender/MPstats and ROI files must be in one directory. Note this
% will take significantly longer to run on MPStats files that have not been
% converted to MPRender. Consider running convert_MPStatstoMPRender.m first
%
% Alex Settle & Miguel de Jesus
% Memorial Sloan Kettering Cancer Center
% Morgan Huse Laboratory, Department of Immunology
% 2023

%% Choose the directory containing all of the data, set options
file_path = uigetdir;
file_list = dir(file_path);

out_path =  uigetdir;
addpath(out_path);


plot_choice = questdlg(['Would you like to make results overview pages' ...
    ' for every calculated synapse?'],'Plot Choice','Yes','No','No');

verbose = 1;

for i = 1:length(file_list)
    if isempty(strfind(file_list(i).name,'MPRender'))
           continue
    end
    disp(['Loading...' file_list(i).name ' ... ' datestr(datetime)])  
    file_name = file_list(i).name;
    index_match = strfind(file_name,'MPRender');
    match_2 = index_match(end);
    roi_filename = [file_name(1:match_2-1),'ROI',file_name(match_2+length('MPRender'):end)];

    MPhandle = load(fullfile(file_path,file_name));
    ROIhandle = load(fullfile(file_path,roi_filename));
    MPStats = MPhandle.MPStats;
    %to smooth any inconsistencies in variable naming
    roidata_Field = fieldnames(ROIhandle);
    roidata = ROIhandle.(roidata_Field{1});

    MP_size = size(MPStats,2);
    ROI_size = size(roidata,2);
    if MP_size ~= ROI_size
        error(['Error: MPStats and ROI do not match in size: '...
            file_name '  //  ' roi_filename])
    end
    synapseResults = initializeSynapseResultsStruct();
    for j = 1:MP_size
        if verbose == 1
            disp(['Calculating ... ' MPStats(j).FileName ' ... ' datestr(datetime)])
        end
        synapseResults(j) = calculateSynapseStats(MPStats(j),roidata(j));
        if strcmp(plot_choice,'Yes')
            outpng = strrep(MPStats(j).FileName,'.tif',' - ResultsOverview.png');
            generateSynapseResultsPage(MPStats(j),roidata(j),synapseResults(j));
            saveas(gcf,fullfile(out_path,outpng))
            close
        end

    end
    outmat = strrep(file_name,'MPRender','Results');
    save(fullfile(out_path,outmat),'synapseResults');
end