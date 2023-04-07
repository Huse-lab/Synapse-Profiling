%% Analyze curvature distance data (V4.2)

% 2021-02-11: Wrote code that overlays circular bins on top of all synapses
% and then takes the data that fall within those bins. Morgan prefers a
% method in which the synapse edge is always 1.0 at any angle. Will try
% hard-coding this.

% 2021-03-30: Added assign_celltype.m to section 1. 2021-04-01...03: (v2.0)
% Streamlined this code using functions. Combined this code with
% analyze_locationstrengthmatrix_fromdistancedata.m and
% analyze_celltypeSpecific_averageradialPDF.m in order to separately
% calculate distance-curvature matrices, then combine them as necessary.

% 2021-04-05: Adjusted the y-axis slightly so they are aligned with their colored data,
% by padding an extra  row on top of the D-C matrices (full of 0s).

% 2021-04-19: Upon Julie's suggestion, re-binned the curvature (y-) axis to
% remove sharp dropoffs, and just to represent the data more faithfully.
% Made adjustments to calculate_DistanceCurvatureMatrix function and plot_DCM to take in
% any arbitrary array of numbers to use as the bins.

% 2021-05-30: Added options for Gaussian Curvature using a different magnitude regime
% for rendering.

% 2021-08-05: Added intensity profiling, and then streamlined the
% code to handle MC, GC, and cell stain data.

%% Initialization: compile distance_data structure into a single structure for analysis

compiled_data = struct('Filename',{},'Celltype',{},'Distance_NormalizedtoEdge',[],'MeanCurvature',[],'GaussianCurvature',[],'CellStain',[]);

path = uigetdir(pwd,'Select folder containing distance data for radial profiling.');
figlist = dir(fullfile(path, '*distance_data*'));

% Loop to fetch data from each data structure and stuff it into
% compiled_data
for idx=1:length(figlist)
    filename = fullfile(path,figlist(idx).name);
    load(filename);
    compiled_data(idx).Filename = figlist(idx).name;
    compiled_data(idx).Distance_NormalizedtoEdge = distance_data.NormalizedRadialD;
    if isfield(distance_data,'MeanCurvature')
        compiled_data(idx).MeanCurvature = distance_data.MeanCurvature;
    end
    if isfield(distance_data,'GaussianCurvature')
        compiled_data(idx).GaussianCurvature = distance_data.GaussianCurvature;
    end    
    if isfield(distance_data,'StainInt')
        compiled_data(idx).CellStain = distance_data.StainInt;
    end
end

 compiled_data = assign_celltype(compiled_data);
 compiled_data = sortStruct(compiled_data,'Celltype');
   
%% Generate distance-curvature matrices

% Find indices of every cell type in the compiled_data structure.
unique_idx = unique_datainStructure(compiled_data);

% Streamline radial data compilation and distance-curvature matrix
% generation
save_choice = 2;
    while ~any(save_choice==[0 1])
        save_choice = input("Save new DC matrix files? 0: No, just calculate. 1: Yes, save new files. ");
    end

nbin = input("How many radial bins? ");
 
% Separately tackle the calculations depending on which signals exist.
if isfield(distance_data,'MeanCurvature')
    cbin = -1:0.1:1;
    input_mat = struct('Distance_NormalizedtoEdge',[],'MeanCurvature',[]);
    for idx = 1:size(unique_idx,1)
        % Use a counter variable to organize the indices better, especially for
        % compiling data of a given type.
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            input_mat.Distance_NormalizedtoEdge = compiled_data(jdx).Distance_NormalizedtoEdge;
            input_mat.MeanCurvature = compiled_data(jdx).MeanCurvature;
            dummy_DCM = calculate_DistanceCurvatureMatrix(input_mat,nbin,cbin);    
            compiled_data(jdx).RadialProfile_MeanCurvature = dummy_DCM;    
            clear dummy_DCM input_mat;
        end
    end  
end

if isfield(distance_data,'GaussianCurvature')
    cbin = -0.3:0.03:0.3;
    % Keep MeanCurvature name for function
    input_mat = struct('Distance_NormalizedtoEdge',[],'MeanCurvature',[]);
    for idx = 1:size(unique_idx,1)
        % Use a counter variable to organize the indices better, especially for
        % compiling data of a given type.
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            input_mat.Distance_NormalizedtoEdge = compiled_data(jdx).Distance_NormalizedtoEdge;
            input_mat.MeanCurvature = compiled_data(jdx).GaussianCurvature;
            dummy_DCM = calculate_DistanceCurvatureMatrix(input_mat,nbin,cbin);    
            compiled_data(jdx).RadialProfile_GaussianCurvature = dummy_DCM;    
            clear dummy_DCM input_mat;
        end
    end  
end

if isfield(distance_data,'StainInt')
    cbin = 0:0.05:1;
    % Keep MeanCurvature name for function
    input_mat = struct('Distance_NormalizedtoEdge',[],'MeanCurvature',[]);
    for idx = 1:size(unique_idx,1)
        % Use a counter variable to organize the indices better, especially for
        % compiling data of a given type.
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            input_mat.Distance_NormalizedtoEdge = compiled_data(jdx).Distance_NormalizedtoEdge;
            % Normalize the pixel intensities to [0,1] range.
            stain_int = compiled_data(jdx).CellStain;
                input_mat.MeanCurvature = rescale(stain_int);
            dummy_DCM = calculate_DistanceCurvatureMatrix(input_mat,nbin,cbin);    
            compiled_data(jdx).RadialProfile_CellStain = dummy_DCM;    
            clear dummy_DCM input_mat;
        end
    end  
end

% Save data
if save_choice==1
    for kdx=1:length(compiled_data)
        matrix_name = extractBefore(compiled_data(kdx).Filename,'_distance');
        if isfield(compiled_data,'RadialProfile_MeanCurvature')
            savename = strcat(matrix_name,'_radprof_MC.mat');
            radialprofile_MC = compiled_data(kdx).RadialProfile_MeanCurvature;
            save(savename,'radialprofile_MC');
        end
        if isfield(compiled_data,'RadialProfile_GaussianCurvature')
            savename = strcat(matrix_name,'_radprof_GC.mat');
            radialprofile_GC = compiled_data(kdx).RadialProfile_GaussianCurvature;
            save(savename,'radialprofile_GC');
        end
        if isfield(compiled_data,'RadialProfile_CellStain')
            savename = strcat(matrix_name,'_radprof_StainInt.mat');
            radialprofile_stainInt = compiled_data(kdx).RadialProfile_CellStain;
            save(savename,'radialprofile_stainInt');
        end
    end
end
 
%% Population-average the radial profiles
% Adapted from analyze_celltypeSpecific_averagedradialPDF.m

if exist('compiled_data','var')
    fprintf("Compile the data by type in order to calculate a population PDF.");
    unique_idx = unique_datainStructure(compiled_data); %Use "Celltype"
end

pop_radprof = struct('Population_subtype',{});

if isfield(compiled_data,'RadialProfile_MeanCurvature')
    for idx = 1:size(unique_idx,1)
        pop_radprof(idx).Population_subtype = unique_idx{idx,1};
        dummy = [];
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            dummy(:,:,counter) = compiled_data(jdx).RadialProfile_MeanCurvature;
        end
        cum_mat = [];
        cum_mat = sum(dummy,3) / size(dummy,3);
        pop_radprof(idx).MC = cum_mat;
    end
end
if isfield(compiled_data,'RadialProfile_GaussianCurvature')
    for idx = 1:size(unique_idx,1)
        pop_radprof(idx).Population_subtype = unique_idx{idx,1};
        dummy = [];
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            dummy(:,:,counter) = compiled_data(jdx).RadialProfile_GaussianCurvature;
        end
        cum_mat = [];
        cum_mat = sum(dummy,3) / size(dummy,3);
        pop_radprof(idx).GC = cum_mat;
    end
end
if isfield(compiled_data,'RadialProfile_CellStain')
    for idx = 1:size(unique_idx,1)
        pop_radprof(idx).Population_subtype = unique_idx{idx,1};
        dummy = [];
        for counter = 1:unique_idx{idx,3}
            jdx = unique_idx{idx,2} + counter - 1;
            dummy(:,:,counter) = compiled_data(jdx).RadialProfile_CellStain;
        end
        cum_mat = [];
        cum_mat = sum(dummy,3) / size(dummy,3);
        pop_radprof(idx).StainInt = cum_mat;
    end
end
        
%% Display all profiles per category

if ~exist('nbin','var')
    nbin = input("How many radial bins? ");
end

plot_type = 3;
while ~any(plot_type==[0 1 2])
    plot_type = input("Which data to plot? 0: Mean Curvature, 1: Gaussian Curvature, 2: Cell Stain ");
end

for idx=1:length(pop_radprof)
    figure
        if plot_type==0
            cbin = -1:0.1:1;
            mat_input = pop_radprof(idx).MC;
        elseif plot_type==1
            cbin = -0.3:0.03:0.3;
            mat_input = pop_radprof(idx).GC;
        elseif plot_type==2
            cbin = 0:0.05:1;
            mat_input = pop_radprof(idx).StainInt;
        end
    plot_DCM(mat_input,pop_radprof(idx).Population_subtype,cbin);
end

%% Pairwise difference heatmap plotting: name any variable for plotting mat_trans

if ~exist('nbin','var')
    nbin = input("How many radial bins? ");
end

plot_type = 3;
while ~any(plot_type==[0 1 2])
    plot_type = input("Which data to plot? 0: Mean Curvature, 1: Gaussian Curvature, 2: Cell Stain ");
end

if plot_type==0
    cbin = -1:0.1:1;
elseif plot_type==1
    cbin = -0.3:0.03:0.3;
elseif plot_type==2
    cbin = 0:0.05:1;
end

figure
plot_diff_DCM(mat_trans,cbin)

%% Functions

% Generate distance curvature matrix from compiled distance_data
% (compiled_data).
% Updated 2021-04-19: REQUIRES the number of radial bins (nbin) and an array of
% curvature bin edges (cbin) as inputs.
function mat_trans = calculate_DistanceCurvatureMatrix(data_input,nbin,cbin)

 % Generate linearized distance and curvature arrays from the
 % cell-type-specific data structures.
        distances = data_input.Distance_NormalizedtoEdge(:)';
        curvatures = data_input.MeanCurvature(:)';
    % Bin the arrays by radial distance. Define max as 1.0 to normalize
    % across cells.
        edges = linspace(0,1,(nbin+1));
        curvatures(2,:) = discretize(distances,edges);
    % Bin by curvatures
        curvature_window = struct('Radial_bin',[],'Curvatures',[],'Curvature_bin_freqs',[]);
        for idx=1:nbin
            curvature_window(idx).Radial_bin = idx;
                c_inShell = curvatures(1,curvatures(2,:)==idx);
                % Calculate frequencies of curvatures per bin/shell
                c_inShell(2,:) = discretize(c_inShell,cbin);
                c_binned = [1:length(cbin)-1];
                c_freq = histc(c_inShell(2,:),c_binned);
            curvature_window(idx).Curvatures = c_inShell(1,:);
            curvature_window(idx).Curvature_bin_freqs = [c_freq];
        end
        clear c_inShell c_binned c_freq
    % Distance-curvature matrix (Location-strength binning)
        c_length = length(curvature_window(1).Curvature_bin_freqs);
        mat=zeros(nbin,c_length);
            for idx = 1:nbin
                for jdx = 1:c_length
                    mat(idx,jdx) = curvature_window(idx).Curvature_bin_freqs(jdx);
                end
            end
        mat_norm=zeros(nbin,c_length);
            for idx=1:nbin
                mat_norm(idx,:)=mat(idx,:)/sum(mat(idx,:));
            end
        mat_norm(isnan(mat_norm(:,:)))=0;
        mat_trans=mat_norm';
        
end

% Plot the DC matrix in my preferred format.
function plot_DCM(input_dcm,input_title,cbin)
 
    nbin = size(input_dcm,1);
    % For padded margins
        input_dcm = [input_dcm;zeros(1,nbin)];
    surf(input_dcm), view(2)
    % For generating y-intervals based on the input array
        cbin_string = string(cbin);
        cbin_tick = cell(1,length(cbin)-1);
            for idx=1:length(cbin_tick)
                tickstring = strcat('[',cbin_string(idx),':',cbin_string(idx+1),']');
                cbin_tick{idx} = tickstring;
            end
        set(gca,'YTickLabel',cbin_tick)  
    % Set all other parameters
        set(gca,'XLim',[1 nbin],'XTick',[1.5 0.5+nbin/2 nbin-0.5],'XTickLabel',{'0' '0.5' '1.0'})       
        set(gca,'YLim',[1 length(cbin)],'YTick',[1.5:1:length(cbin)+0.5])
        xlabel('Radial distance')
        % Guess parameter for labeling based on cbin limits
        if min(cbin) == -1 && max(cbin) == 1
            parameter_label = 'Local mean curvature';
            caxis([0 0.4])
        elseif min(cbin) == -0.3 && max(cbin) == 0.3
            parameter_label = 'Local Gaussian curvature';
            caxis([0 0.1])
        elseif min(cbin) == 0 && max(cbin) == 1
            parameter_label = 'Cell stain intensity';
        end
        titlestring = strcat("Radial profile: ",input_title);
        title(titlestring)
        colorcet('L17')
        cb=colorbar();
        cb.Title.String='pdf';
        shading flat
        set(gcf,'color',[1 1 1])
        daspect([1 1 1])
        ylabel(parameter_label)
    
end

% Plot the DC matrix in my preferred format.
function plot_diff_DCM(input_dcm,cbin)
    
    string_celltypes = input("What cells are we comparing? ",'s');
    nbin = size(input_dcm,1);
    % For padded margins
        input_dcm = [input_dcm;zeros(1,nbin)];
    surf(input_dcm), view(2)
    % For generating y-intervals based on the input array
        cbin_string = string(cbin);
        cbin_tick = cell(1,length(cbin)-1);
            for idx=1:length(cbin_tick)
                tickstring = strcat('[',cbin_string(idx),':',cbin_string(idx+1),']');
                cbin_tick{idx} = tickstring;
            end
        set(gca,'YTickLabel',cbin_tick)  
    % Set all other parameters        
        set(gca,'XLim',[1 nbin],'XTick',[1.5 0.5+nbin/2 nbin-0.5],'XTickLabel',{'0' '0.5' '1.0'})       
        set(gca,'YLim',[1 length(cbin)],'YTick',[1.5:1:length(cbin)+0.5])
        xlabel('Radial distance')
        if min(cbin) == -1 && max(cbin) == 1
            ylabel('Local mean curvature');
        elseif min(cbin) == -0.3 && max(cbin) == 0.3
            ylabel('Local Gaussian curvature');
        elseif min(cbin) == 0 && max(cbin) == 1
            ylabel('Cell stain intensity');
        end
        titlestring = strcat("Differential radial profile: ",string_celltypes);
        title(titlestring);
        colorcet('D9')
        cb=colorbar();
        cb.Title.String='\Delta pdf';
        shading flat
        set(gcf,'color',[1 1 1])
        daspect([1 1 1])
        caxis([-0.05 0.05])
    
end
