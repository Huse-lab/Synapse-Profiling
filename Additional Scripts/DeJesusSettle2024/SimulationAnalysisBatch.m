%% Initialization

file_path = uigetdir(pwd,"Where is the .csv data stored? ");
file_list = dir(fullfile(file_path,'*.csv'));

%% Loop through all csvs, calculate curvature profiles and stats

T = table();
k=1;

for i = 1:length(file_list)
    i
    file = file_list(i);
    if startsWith(file.name,'z_')
        tic
        file.name
        k
        param_parse = strsplit(strrep(file.name,'z_','param_'),'_rep_');
        param_file = [param_parse{1} '.csv'];
        act_file = strrep(file.name,'z_','prot_');
        metadata = struct2table(read_parameter_file(param_file),'AsArray',true);
        z = csvread(file.name);
        act_data = csvread(act_file);
        curvatureData = struct2table(calculateSimulatedCurvatureStats(z,act_data),'AsArray',true);
        table_row = [metadata curvatureData];
        T = [T;table_row];
        k=k+1;
        toc
    end
end

%% clean up table to make saving simpler for plotting.

cleanT = removevars(T,{'concave_area_map','DegranMap','mc_all','convex_area_map'});
cleanT = sortrows(cleanT,{'NumProtrusions','MeanClusterDiameter'});
writetable(cleanT,'05_05_23_ResultsTable.csv')


%% Average replicates

%metadata labels
meta_fields = {'ClusterSize','NumProtrusions','SurfaceTension','BendingRigidity'};
metafields_wname = [meta_fields 'FileName'];

unique_field_array = {};
for i = 1:length(meta_fields)
    unique_field_array = [unique_field_array; unique(cleanT.(meta_fields{i}))'];
end
unique_combos = combvec(unique_field_array{:});

sumTable = table();
warning('off','all');
for i = 1:size(unique_combos,2)
    for j = 1:length(meta_fields)
        sumTable.(meta_fields{j})(i) = unique_combos(j,i);
    end
end

results_fields = setdiff(cleanT.Properties.VariableNames,metafields_wname);

for k = 1:length(results_fields)
    field_name = [results_fields{k} ' Mean'];
    field_std = [results_fields{k} ' Stdv'];
    mean_values = zeros(height(sumTable),1);
    std_devs = zeros(height(sumTable),1);
    for j = 1:height(sumTable)
        replicates = cleanT((sum(cleanT{:,meta_fields} == unique_combos(:,j)',2))==length(meta_fields),:);
        mean_values(j) = mean(replicates.(results_fields{k}));
        std_devs(j) = std(replicates.(results_fields{k})); 
    end
    newColumns = table(mean_values,std_devs,'VariableNames',{field_name,field_std});
    sumTable = [sumTable newColumns];
end

writetable(sumTable,'Results_AveragedTable.csv')
