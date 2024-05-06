%Script for running Zernicke rotational transformations on 50x50 gridded
%synapses, cleaned up for easier use by AHS 9/15/21

% Edited by MDJ 2021-11-30 for compatibility with new data formats/names.
% Copied over _cleanscript_v2 from Alex Code Drop on 2021-12-14
% Combined with Alex/Zernike_transform_cleanscript_v3.m on 2022-01-10

% Added from _cleanscript_v4 from Alex Code Drop on 2022-06-15: fixed array
% recording to include both Rotationally collapsed and un-collapsed
% Z-arrays

% 2022-10-28 Added functions to better identify synapse indices within
% certain regions by polygon drawing

% 2022-10-30: V4, added k-means functions to split the Z-PCA
% 2023-12-22: V4(again), restructured the code to handle multiple types of
% grids in parallel. Important for plotting MC/actin plots side-by-side in
% Z-PCA-gating

%% Initialization

Grid_type_mc = 1;
Grid_type_gc = 0;
Grid_type_int = 0;

Zernike_order = 15;

%% Set the color scheme for subsequent plots here
    % kulay = brewermap(length(unique(celltype_vect)),'set1');
% 10-color for simulations
kulay = brewermap(10,'paired');
% 2-color Blue-Red for BMDM vs. OT-1
kulay = brewermap(2,'set1');
% 8-color for the simulations
kulay = brewermap(8,'dark2');
% 16-color for simulations at Surface Tension = 1 & 2
kulay = [brewermap(8,'set2');brewermap(8,'dark2')];
% 4-color for Tn/Teff/Tmem/Texh experiment
kulay = brewermap(8,'dark2');
kulay = kulay(3:6,:);
% 4-color for sgNT/Talin/WASp/WAVE2 + non-synapse: 
kulay = brewermap(8,'set1');
kulay = kulay([3:5,7],:); kulay = [kulay;0 0 0];
% 3-color for 300/500/1000
kulay = [128, 128, 128;0, 0, 255;249 64 64]/255;
% 3-color for sgNT/PTEN/WASp
kulay = [128, 128, 128;85, 160, 251;255, 128, 128]/255;

% 4-color for Tn/Teff/Tmem/Texh experiment
kulay = brewermap(8,'dark2');
kulay = kulay(3:6,:);

% 2-color for sgNT/sgPTEN at 300Pa
kulay = [0.4 0.4 0.4;1 0 0];
% try x's and o's for the stiffness

% 2-color Gray-Pink for CD4 vs. CD8:
kulay = [153 153 153;247, 129, 191]/255;

%% Import grids into a cell array 

% Pick folders with grid data
if Grid_type_mc
    mc_file_path = uigetdir(pwd,'Please select folder with MC grids');
    addpath(mc_file_path);
end
if Grid_type_gc
    gc_file_path = uigetdir(pwd,'Please select folder with GC grids');
    addpath(gc_file_path);
end
if Grid_type_int
    int_file_path = uigetdir(pwd,'Please select folder with cell stain grids');
    addpath(int_file_path);
end

%% Calculate Zernike coefficients

% Initialize dummy arrays
mc_synapses = {};gc_synapses = {};int_synapses = {};
names_mc = {};names_gc = {};names_int = {};

% Calculate the Zernike coefficients and store them into data structures
if Grid_type_mc
    % Load grids
    files = dir(mc_file_path);
    for file = files'
        if contains(file.name,['_MeanC'])
            load(file.name)
            mc_synapses = [mc_synapses; mc_i];
            names_mc = [names_mc;file.name];
        end
    end
    % Calculate Zernike coefficients
    [C_array_mc,a_array_mc,Mr_mc,Nr_mc] = zernikulator(mc_synapses,0,Zernike_order);
    % Generate data structure and store data
    coeff_struct_mc = struct('FileName',{},'C_zernike',{},'C_adjusted_zernike',{});
    for k = 1:length(names_mc)
        coeff_struct_mc(k).FileName = names_mc(k);
        coeff_struct_mc(k).C_zernike = a_array_mc(k,:);
        coeff_struct_mc(k).C_adjusted_zernike = C_array_mc(k,:);
    end
    coeff_struct_mc = assign_celltype(coeff_struct_mc);
    coeff_struct_mc = sortStruct(coeff_struct_mc,'Celltype');
end
if Grid_type_gc
    files = dir(gc_file_path);
    for file = files'
        if contains(file.name,['_GaussC'])
            load(file.name)
            mc_synapses = [gc_synapses; gc_i];
            names_gc = [names_gc;file.name];
        end
    end
    [C_array_gc,a_array_gc,Mr_gc,Nr_gc] = zernikulator(gc_synapses,0,Zernike_order);
    coeff_struct_gc = struct('FileName',{},'C_zernike',{},'C_adjusted_zernike',{});
    for k = 1:length(names_gc)
        coeff_struct_gc(k).FileName = names_gc(k);
        coeff_struct_gc(k).C_zernike = a_array_gc(k,:);
        coeff_struct_gc(k).C_adjusted_zernike = C_array_gc(k,:);
    end
    coeff_struct_gc = assign_celltype(coeff_struct_gc);
    coeff_struct_gc = sortStruct(coeff_struct_gc,'Celltype');
end
if Grid_type_int
    files = dir(int_file_path);
    for file = files'
        if contains(file.name,['_Cellstain'])
            load(file.name)
            % Can normalize per file or globally. If globally, remove the next
            % line and then divide int_synapses by its max at the end of
            % the forloop.
                 int_i = int_i ./ max(max(int_i));
                int_synapses = [int_synapses; int_i];
            names_int = [names_int;file.name];
        end
    end
    [C_array_int,a_array_int,Mr_int,Nr_int] = zernikulator(int_synapses,0,Zernike_order);
    coeff_struct_int = struct('FileName',{},'C_zernike',{},'C_adjusted_zernike',{});
    for k = 1:length(names_int)
        coeff_struct_int(k).FileName = names_int(k);
        coeff_struct_int(k).C_zernike = a_array_int(k,:);
        coeff_struct_int(k).C_adjusted_zernike = C_array_int(k,:);
    end
    coeff_struct_int = assign_celltype(coeff_struct_int);
    coeff_struct_int = sortStruct(coeff_struct_int,'Celltype');
end

%% Manual: manually annotate celltypes if needed, sortStruct, and save

%% Calculate & plot Z-PCA

list_of_datatypes = {...
    'MC',...                % Choice 1
    'GC',...                % Choice 2
    'Cellstain',...         % Choice 3
    'MC + Cellstain',...    % Choice 4
    'GC + Cellstain',...    % Choice 5
    'MC + GC',...           % Choice 6
    'all'...                % Choice 7
    };
[choice_vector,~] = listdlg('Promptstring',{'What data should Z-PCA sort on?'},'ListString',list_of_datatypes,'SelectionMode','single');
warning('off','all');

feature_vect = [];
celltype_vect = {};

% Depending on the choice made, Z-PCA will be calculated based on 1 or a
% combination of the gridded data types. Will raise error if there is a
% mismatch with the data types available.
switch choice_vector
    case 1      % MC only
        for i = 1:length(coeff_struct_mc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_mc(i).Celltype)];
            feature_vect = [feature_vect; coeff_struct_mc(i).C_adjusted_zernike];
        end
    case 2      % GC only
        for i = 1:length(coeff_struct_gc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_gc(i).Celltype)];
            feature_vect = [feature_vect; coeff_struct_gc(i).C_adjusted_zernike];
        end
    case 3      % Cellstain only
        for i = 1:length(coeff_struct_int)
            celltype_vect = [celltype_vect; num2str(coeff_struct_int(i).Celltype)];
            feature_vect = [feature_vect; coeff_struct_int(i).C_adjusted_zernike];
        end
    case 4      % MC + Cellstain
        for i = 1:length(coeff_struct_mc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_mc(i).Celltype)];
            feature_vect = [feature_vect; [coeff_struct_mc(i).C_adjusted_zernike coeff_struct_int(i).C_adjusted_zernike]];
        end
    case 5      % GC + Cellstain
        for i = 1:length(coeff_struct_gc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_gc(i).Celltype)];
            feature_vect = [feature_vect; [coeff_struct_gc(i).C_adjusted_zernike coeff_struct_int(i).C_adjusted_zernike]];
        end
    case 6      % MC + GC
        for i = 1:length(coeff_struct_mc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_mc(i).Celltype)];
            feature_vect = [feature_vect; [coeff_struct_mc(i).C_adjusted_zernike coeff_struct_gc(i).C_adjusted_zernike]];
        end
    case 7      % All
        for i = 1:length(coeff_struct_mc)
            celltype_vect = [celltype_vect; num2str(coeff_struct_mc(i).Celltype)];
            feature_vect = [feature_vect; [coeff_struct_mc(i).C_adjusted_zernike coeff_struct_gc(i).C_adjusted_zernike coeff_struct_int(i).C_adjusted_zernike]];
        end
end

figure
fig_idx_ZPCA = get(gcf,'Number'); % Get the record of this figure so that the Z-PCA gallery can always find it.
hold on
[coeff, score, latent, tsquared, explained] = pca(feature_vect);
gscatter(score(:,1),score(:,2),celltype_vect,kulay,'.',20)
daspect([1 1 1])
xlabel(['PC 1 (' num2str(explained(1)) '%)'])
ylabel(['PC 2 (' num2str(explained(2)) '%)'])
set(gca,'TickLength',[0 0])
box off

%% Manually define regions of the Z-PCA plot

% Start with the Z-PCA figure

gates_ZPCA = struct('ROI_info',{},'Names_contained',{});
choice=1;

while choice==1
    
    draw_decision = questdlg('Would you like to draw a new polygon?','Yes','No');
    switch draw_decision
        case 'Yes'
        % If the choice is "Yes", keep going!
        figure(fig_idx_ZPCA)
           % Draw a polygon to define an ROI
            roi_ZPCA = drawpolygon('FaceAlpha',0,'Color',[0 0 0],'StripeColor',[1 1 1],'LineWidth',1);
            PC_array = [score(:,1) score(:,2)];
            % Accounts for re-sorting by going back to pre-sorted array "names" (AHS)
            within_idx = inpolygon(PC_array(:,1),PC_array(:,2),roi_ZPCA.Position(:,1),roi_ZPCA.Position(:,2));
            % Plot galleries by grid type
            if Grid_type_mc
                % Find the names & grids that correspond to the items within the polygon.
                within_names = [coeff_struct_mc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_mc,within_names);
                display_synapses = mc_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                    caxis([-0.5 0.5]);
                    colormap(turquoisebrown);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    box off;
                    axis off;
                end
            end
            if Grid_type_gc
                within_names = [coeff_struct_gc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_gc,within_names);
                display_synapses = gc_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                    caxis([-0.3 0.3]);
                    colormap(greenmagenta);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    title(within_names(idx))
                    box off;
                    axis off;
                end
            end
            if Grid_type_int
                within_idx = inpolygon(PC_array(:,1),PC_array(:,2),roi_ZPCA.Position(:,1),roi_ZPCA.Position(:,2));
                within_names = [coeff_struct_int(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_int,within_names);
                display_synapses = int_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                        stainmap = ones(256,3)*.05;
                        stainmap(:,2) = linspace(.05,.9,256);
                        caxis([0 max(max(display_synapses{idx}))]);
                        colormap(stainmap);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    box off;
                    axis off;
                end
            end
            % Store info into the structure
            gates_ZPCA(length(gates_ZPCA)+1).ROI_info = roi_ZPCA;
            gates_ZPCA(length(gates_ZPCA)).Names_contained = within_names;
            
        case 'No'
            break
    end
        
end

for idx=1:length(gates_ZPCA)
    gates_ZPCA(idx).Polygon_coords = gates_ZPCA(idx).ROI_info.Position;
end

%% Replot manual Z-PCA gates with clean polyshapes 

% If replotting later on, you can use a for-loop through the roi_ZPCA
% contents and then plot via

for idx=1:length(gates_ZPCA)
 patch(gates_ZPCA(idx).Polygon_coords(:,1),gates_ZPCA(idx).Polygon_coords(:,2),[1 1 1],'FaceAlpha',0,'LineStyle','--');
end

%% Automatically split Z-PCA plot into quantiles

% Define the quantiles for grouping patterns into n-groups by PC1 value. 
% Follows n+1 pattern: e.g., for quartiles, n_quantiles = 3.
n_quantiles = 2;

gates_ZPCA = struct('ROI_info',{},'Names_contained',{});
choice=1;
fig_idx_ZPCA = get(gcf,'Number');
if n_quantiles > 1
    split_x = quantile(score(:,1),n_quantiles);
elseif n_quantiles==1
    split_x = median(score(:,1));
end
bounds_x = get(gca,'xlim');
regions_x = [bounds_x(1) split_x bounds_x(2)];

for idx=1:n_quantiles+1
    if idx<n_quantiles+1
        xline(split_x(idx),'-','Color',[0 0 1],'LineWidth',2.5);
    end
end

% Ask user to plot gallery or not:
list = {'Yes, plot the galleries','No, aint nobody got time for that'};
[choice_plot,~] = listdlg('Promptstring',{'Do you want to plot individual Region galleries? (Takes time)'},'ListString',list,'SelectionMode','single');
warning('off','all');

if choice_plot==1
    for idx=1:n_quantiles+1
        within_idx = score(:,1) >= regions_x(idx) & regions_x(idx+1) > score(:,1);
        switch choice_vector
            case 1
                within_names = [coeff_struct_mc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_mc,within_names);
                display_synapses = mc_synapses(display_idx);
            case 2
                within_names = [coeff_struct_gc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_gc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_gc,within_names);
                display_synapses = gc_synapses(display_idx);
            case 3
                within_names = [coeff_struct_int(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_int);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_int,within_names);
                display_synapses = int_synapses(display_idx);
        end
        % Store info into the structure
        gates_ZPCA(length(gates_ZPCA)+1).ROI_info = [regions_x(idx) regions_x(idx+1)];
        gates_ZPCA(length(gates_ZPCA)).Names_contained = within_names;            

            % Plot galleries by grid type
            if Grid_type_mc
                % Find the names & grids that correspond to the items within the polygon.
                within_names = [coeff_struct_mc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_mc,within_names);
                display_synapses = mc_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    if contains(within_names(idx),'5CC7')
                        title('4')
                    elseif contains(within_names(idx),'OT')
                        title('8')
                    end
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                    caxis([-0.5 0.5]);
                    colormap(turquoisebrown);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    box off;
                    axis off;
                end
            end
            if Grid_type_gc
                within_names = [coeff_struct_gc(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_gc,within_names);
                display_synapses = gc_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                    caxis([-0.3 0.3]);
                    colormap(greenmagenta);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    title(within_names(idx))
                    box off;
                    axis off;
                end
            end
            if Grid_type_int
                within_names = [coeff_struct_int(within_idx).FileName];
                % [~,display_idx] = ismember(within_names,names_mc);
                % ^changed line to below due to update of ismember function
                display_idx = ismember(names_int,within_names);
                display_synapses = int_synapses(display_idx);
                % Create gallery of the synapses collected in ROI
                n_grid = length(display_synapses);
                n_row = ceil(sqrt(n_grid));
                n_col = ceil(n_grid/n_row);
                figure
                gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
                for idx=1:length(display_synapses)
                    axes(gallery(idx));
                    imagesc(display_synapses{idx});
                    view(2);
                    daspect([1 1 1]); pbaspect([1 1 1]);
                    shading flat;
                    set(gca,'Xtick',[]); set(gca,'Ytick',[])
                        stainmap = ones(256,3)*.05;
                        stainmap(:,2) = linspace(.05,.9,256);
                        caxis([0 max(max(display_synapses{idx}))]);
                        colormap(stainmap);
                    set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                    box off;
                    axis off;
                end
            end
%             
%     % Create gallery of the synapses collected in ROI
%         n_grid = length(display_synapses);
%         n_row = ceil(sqrt(n_grid));
%         n_col = ceil(n_grid/n_row);
%         figure
%         gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
% 
%         for jdx=1:length(display_synapses)
%                 axes(gallery(jdx));
%                 surf(display_synapses{jdx});
%                 view(2);
%                 daspect([1 1 1]); pbaspect([1 1 1]);
%                 shading flat;
%                 set(gca,'Xtick',[]); set(gca,'Ytick',[])
%                 colormap(turquoisebrown);
%                 caxis([-0.3 0.3])
%                 set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
%                 % titlestring = strcat(griddata_struct(idx).Celltype," (n = ",num2str(n_grid),")");
%                 % sgtitle(titlestring);
%         end

        for jdx=1:n_row*n_col
            axes(gallery(jdx));
            axis off;
        end

    end
else
    for idx=1:n_quantiles+1
        within_idx = score(:,1) >= regions_x(idx) & regions_x(idx+1) > score(:,1);
        switch choice_vector
            case 1
                within_names = [coeff_struct_mc(within_idx).FileName];
            case 2
                within_names = [coeff_struct_gc(within_idx).FileName];
            case 3
                within_names = [coeff_struct_int(within_idx).FileName];
        end
        % Store info into the structure
        gates_ZPCA(length(gates_ZPCA)+1).ROI_info = [regions_x(idx) regions_x(idx+1)];
        gates_ZPCA(length(gates_ZPCA)).Names_contained = within_names;            
    end
end

%% Automatically split Z-PCA by k-means

n_kbins = 4;
idx2Region = kmeans(score(:,[1 2]),n_kbins);
% Follows naming order of coeff_struct when Z-PCA was calculated.

% Ask user to plot gallery or not:
list = {'Yes, plot the galleries','No, aint nobody got time for that'};
[choice_plot,~] = listdlg('Promptstring',{'Do you want to plot individual Region galleries? (Takes time)'},'ListString',list,'SelectionMode','single');
warning('off','all');

gates_ZPCA = struct('Polygon_coords',{},'Names_contained',{});

for bin = 1:n_kbins
    within_bin = idx2Region==bin;
    if Grid_type_mc
        synapses_in_bin = coeff_struct_mc(within_bin);
        gates_ZPCA(bin).Names_contained = [synapses_in_bin.FileName];
        display_synapses = mc_synapses(within_bin);
        if choice_plot==1
            % Create gallery of the synapses collected in ROI
            n_grid = length(display_synapses);
            n_row = ceil(sqrt(n_grid));
            n_col = ceil(n_grid/n_row);
            figure
            gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
            for idx=1:length(display_synapses)
                axes(gallery(idx));
                imagesc(display_synapses{idx});
                view(2);
                daspect([1 1 1]); pbaspect([1 1 1]);
                shading flat;
                set(gca,'Xtick',[]); set(gca,'Ytick',[])
                caxis([-0.5 0.5]);
                colormap(turquoisebrown);
                set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                box off;
                axis off;
            end
        end
    end
    if Grid_type_int
        synapses_in_bin = coeff_struct_int(within_bin);
        gates_ZPCA(bin).Names_contained = [synapses_in_bin.FileName];
        display_synapses = int_synapses(within_bin);
        if choice_plot==1
            % Create gallery of the synapses collected in ROI
            n_grid = length(display_synapses);
            n_row = ceil(sqrt(n_grid));
            n_col = ceil(n_grid/n_row);
            figure
            gallery = tight_subplot(n_row,n_col,[.01 .005],[.1 .01],[.01 .01]);
            for idx=1:length(display_synapses)
                axes(gallery(idx));
                imagesc(display_synapses{idx});
                view(2);
                daspect([1 1 1]); pbaspect([1 1 1]);
                shading flat;
                set(gca,'Xtick',[]); set(gca,'Ytick',[])
                    stainmap = ones(256,3)*.05;
                    stainmap(:,2) = linspace(.05,.9,256);
                    caxis([0 max(max(display_synapses{idx}))]);
                    colormap(stainmap);
                set(gca,'XLim',[1 50]); set(gca,'YLim',[1 50]);
                box off;
                axis off;
            end
        end
    end
    
    coords_within = [score(within_bin,1) score(within_bin,2)];
    convhull_bin = convhull(coords_within);
    roi_within = [coords_within(convhull_bin,1) coords_within(convhull_bin,2)];
    gates_ZPCA(bin).Polygon_coords = roi_within;
    
end
            
figure, hold on;
plot(score(:,1),score(:,2),'.','MarkerSize',15);
gscatter(score(:,1),score(:,2),idx2Region,brewermap(n_kbins,'Set2'),'.',20)
daspect([1 1 1])
xlabel(['PC 1'])
ylabel(['PC 2'])
set(gca,'TickLength',[0 0])
box off

%% Pattern category disproportionation

switch choice_vector
    case 1
        ref_filename = [coeff_struct_mc.FileName];
        ref_celltype = {coeff_struct_mc.Celltype};
    case 2
        ref_filename = [coeff_struct_gc.FileName];
        ref_celltype = {coeff_struct_gc.Celltype}; 
    case 3
        ref_filename = [coeff_struct_int.FileName];
        ref_celltype = {coeff_struct_int.Celltype}; 
end

unique_celltypes = unique(ref_celltype);
T_regioncomposition = array2table(rand(n_quantiles+1,length(unique_celltypes)),'VariableNames',unique_celltypes);

for idx=1:length(gates_ZPCA)
    
    names_gated = gates_ZPCA(idx).Names_contained;
    [~,within_names] = ismember(names_gated,ref_filename);
    within_celltypes = ref_celltype(within_names);  
    dataintoTable = zeros(1,length(unique_celltypes));
    for jdx=1:length(unique_celltypes)
        count = ismember(within_celltypes,unique_celltypes(jdx));
        dataintoTable(jdx) = sum(count);
    end
    T_regioncomposition(idx,:) = array2table(dataintoTable);
end

T_regioncomposition

%% Highlight cell types separately: scatter

figure('Position',[10 10 1200 1200])
celltype_list = unique(celltype_vect);
% n_row = ceil(sqrt(length(celltype_list)));
% n_col = ceil(length(celltype_list)/n_row);
n_row = 2;
n_col = 2;

for idx=1:length(unique(celltype_vect))
    
    subplot(n_row,n_col,idx)
    cell_labels = celltype_vect;
    cell_labels(~strcmp(celltype_vect,celltype_list{idx})) = {'Other'};
    if idx==1
        colorvect = [0 0 1;0.75 0.75 0.75];
        sz = [12,6];
    else
        colorvect = [0.75 0.75 0.75;0 0 1];
        sz = [6,12];
    end
    h = gscatter(score(:,1),score(:,2),cell_labels,colorvect,'.',sz);
     if idx==1
         legend(h(1),h(1).DisplayName);set(gca,'Children',flipud(get(gca,'Children')))
     else
         legend(h(2),h(2).DisplayName);
     end
    title(celltype_list{idx})
    legend('off');
    
%     xlabel(['PC 1 (' num2str(explained(1)) '%)'])   
%     ylabel(['PC 2 (' num2str(explained(2)) '%)'])
end

sgtitle('Zernike mode PCA (by synapse type)')

%% Highlight cell types separately: contour

figure('Position',[10 10 1200 1200])
celltype_list = unique(celltype_vect);
n_row = ceil(sqrt(length(celltype_list)));
n_col = ceil(length(celltype_list)/n_row);

for idx=1:length(unique(celltype_vect))
    
    subplot(n_row,n_col,idx)
    cell_labels = celltype_vect;
    cell_labels(~strcmp(celltype_vect,celltype_list{idx})) = {'Other'};

    h = gscatter(score(:,1),score(:,2),cell_labels,[0.75 0.75 0.75],'.',6);
   
    hold on
    x = score(strcmp(celltype_vect,celltype_list{idx}),1);
    y = score(strcmp(celltype_vect,celltype_list{idx}),2);
    x_edge = linspace(min(x),max(x),7);
    y_edge = linspace(min(y),max(y),7);
    N_gram = histcounts2(x,y,x_edge,y_edge,'Normalization','countdensity');
    [X,Y] = meshgrid(mean([x_edge(1:end-1);x_edge(2:end)],1),mean([y_edge(1:end-1);y_edge(2:end)],1));
    g = contour(X,Y,N_gram',7);
    % Uncomment out the following indented section of you want the contours
    % forcibly closed. Note this is artifactual but better for
    % visualization.
        x_int = x_edge(2) - x_edge(1); y_int = y_edge(2) - y_edge(1);
        x_edge_pad = [x_edge(1)-x_int/50,x_edge,x_edge(end)+x_int/50];
        y_edge_pad = [y_edge(1)-y_int/50,y_edge,y_edge(end)+y_int/50];
        [X_pad,Y_pad] = meshgrid(mean([x_edge_pad(1:end-1);x_edge_pad(2:end)],1),mean([y_edge_pad(1:end-1);y_edge_pad(2:end)],1));
        N_pad = padarray(N_gram,[1 1],0,'both');
        g = contour(X_pad,Y_pad,N_pad',7);
    colorcet('L8')
    title(celltype_list{idx})
    legend('off');
    xlabel(['PC 1 (' num2str(explained(1)) '%)'])   
    ylabel(['PC 2 (' num2str(explained(2)) '%)'])
    
end
sgtitle('Zernike mode PCA (by synapse type)')

%% Zernike mode breakdown - only works if you run the Zernike calcs on old version Zernike_decomposition_v2.m
% Adapted from AHS's pca_important_modes.m (2022-06-15)

N = []; M = [];
for n = 0:Zernike_order
    N = [N n*ones(1,n+1)];
    M = [M -n:2:n];
end


%initialize circular synapse parameters for zernike calculation 
L = size(mc_synapses{1},1);
X = -1:2/(L-1):1;
[x,y] = meshgrid(X);
x = x(:); y = y(:);
[theta,r] = cart2pol(x,y); 
is_in_circle = ( r <= 1 );
Z = zernfun(N,M,r(is_in_circle),theta(is_in_circle));

    
figure

% Adapted to do PC1_absvalMIN plots (6) and PC2_absvalMAX plots (6)
subplot(5,6,[1:3,7:9,13:15])
bar(coeff(:,1),'EdgeColor','none','FaceColor',[.5 .5 .5]), pbaspect([1 1 1])
ylim([-.5 .5])
title(['PC1 top 6 (' num2str(explained(1)) ' %)'])
subplot(5,6,[4:6,10:12,16:18])
bar(coeff(:,2),'EdgeColor','none','FaceColor',[.5 .5 .5]), pbaspect([1 1 1])
ylim([-.5 .5])
title(['PC2 top 6 (' num2str(explained(2)) ' %)'])

% Find the top ranked polynomials in Z-adj space
[~,top_PC1] = maxk(abs(coeff(:,1)),6); [~,top_PC2] = maxk(abs(coeff(:,2)),6);
topN_1_Zadj = Nr_int(top_PC1); topM_1_Zadj = Mr_int(top_PC1);
topN_2_Zadj = Nr_int(top_PC2); topM_2_Zadj = Mr_int(top_PC2);

count=0;
for i=[19:21,25:27]
    count=count+1;
    Z_orig_idx = intersect(find(N==topN_1_Zadj(count)),find(M==topM_1_Zadj(count)));
    % Alex's method to plot: create an array rec_a that is all =0 except
    % for the desired pattern =1 (at Z_orig_idx), then fetch the pattern
    % from the master Zernike pattern array Z and stuff it into
    % "reconstructed" with the circular mask over it.
    rec_a = zeros(length(a),1);
    rec_a(Z_orig_idx) = 1;
    reconstructed = NaN(size(mc_i));
    reconstructed(is_in_circle) = Z*rec_a;
    subplot(5,6,i)
    imagesc(reconstructed), axis off, pbaspect([1 1 1])
    title([num2str(N(Z_orig_idx)),',',num2str(M(Z_orig_idx))]) 
end

count=0;
for i=[22:24,28:30]
    count=count+1;
    Z_orig_idx = intersect(find(N==topN_2_Zadj(count)),find(M==topM_2_Zadj(count)));
    rec_a = zeros(length(a),1);
    rec_a(Z_orig_idx) = 1;
    reconstructed = NaN(size(mc_i));
    reconstructed(is_in_circle) = Z*rec_a;
    subplot(5,6,i)
    imagesc(reconstructed), axis off, pbaspect([1 1 1])
    title([num2str(N(Z_orig_idx)),',',num2str(M(Z_orig_idx))]) 
end

colorcet('L17')
% text(50,60,'generated by Zernike_decomposition_v2.m',...
%     'Interpreter','None','HorizontalAlignment','Right')

%% Functions

function [C_array,a_array,Mr,Nr] = zernikulator(synapse_array, norm_zern, Zernike_order)

    % Initialize circular synapse parameters for Zernike calculation 
    L = size(synapse_array{1},1);
    X = -1:2/(L-1):1;
    [x,y] = meshgrid(X);
    x = x(:); y = y(:);
    [theta,r] = cart2pol(x,y); 
    N = []; M = [];
    for n = 0:Zernike_order
        N = [N n*ones(1,n+1)];
        M = [M -n:2:n];
    end
    is_in_circle = ( r <= 1 );
    Z = zernfun(N,M,r(is_in_circle),theta(is_in_circle));
    
    % C_array are rotationally adjusted; a_array are specific to the grid
    % orientation (noumenal Zernike coefficients)
    C_array = [];a_array = [];
    % Zernike calculation loop
    for i = 1:length(synapse_array)
        grid_i = synapse_array{i};
        a = Z\grid_i(is_in_circle);
        a_array = [a_array; a'];
        % disp(i);disp(names_mc(i))
        C = [];
        Mr = [];
        Nr = [];
        for j = 1:length(M)
            if M(j) > 0
                c1 = a(j);
                c2 = a(j-M(j));
                Cnew = sqrt(c1^2 + c2^2);
                C = [C Cnew];
                Mr = [Mr M(j)];
                Nr = [Nr N(j)];
            elseif M(j) == 0
                Cnew = a(j);
                C = [C Cnew];
                Mr = [Mr M(j)];
                Nr = [Nr N(j)];
            end
        end
        % Choice to use normalized Zernike coefficients or not (smooth out
        % intensities)
        if norm_zern==1
            C_array = [C_array; normalize(C)];
        elseif norm_zern==0
            C_array = [C_array;C];
        end
    end

end